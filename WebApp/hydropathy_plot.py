from __future__ import annotations
import streamlit as st
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, AutoMinorLocator, FormatStrFormatter
import requests
import numpy as np


######################################################################################## BACKEND ############################################################################


@st.cache()
def Fetch_Sequence_NCBI(Accession_ID: str) -> str:
    """
    :param Accession_ID: Accession ID of the desired protein
    :return: Protein sequence in FASTA format
    """

    url = f"https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?id={Accession_ID}&db=protein&report=fasta&extrafeat=null&conwithfeat=on&hide-cdd=on&retmode=html&withmarkup=on&tool=portal&log$=seqview&maxdownloadsize=1000000"
    Sequence = requests.get(url)
    if Sequence.text.startswith("Failed"):
        return -1
    else:
        return Sequence.text


@st.cache()
def FASTA_Parser(FASTA_seq: str) -> int | str:
    """
    :param FASTA_seq: Protein sequence in FASTA format
    :return: linear aminoacid sequence
    """
    if not FASTA_seq.startswith(">"):
        return -1
    else:
        FASTA_seq = FASTA_seq.split("\n")
        if len(FASTA_seq) < 2:
            return -1
        else:
            AA_seq = "".join(FASTA_seq[1:])
            return AA_seq


@st.cache()
def Preview_Sequence(Seq: str) -> str:
    output = []
    for index in range(0, len(Seq), 10):
        if index % 60 != 0:
            output.append(Seq[index:index + 10])
        else:
            output.append("\n" + Seq[index:index + 10])
    return "    ".join(output)


def Hydropathicity_value_calc(AA_seq_segment: str, Edge_weight=100, Model="Linear Variation") -> float:
    """
    :param AA_seq: Amino acid sequence
    :return: Hydropathy index of the corresponding sequence
    """
    Kyte_Doolittle_scale = {
        "G": -0.4,
        "A": +1.8,
        "V": +4.2,
        "L": +3.8,
        "I": +4.5,
        "P": -1.6,
        "C": +2.5,
        "M": +1.9,
        "S": -0.8,
        "T": -0.7,
        "N": -3.5,
        "Q": -3.5,
        "F": +2.8,
        "Y": -1.3,
        "W": -0.9,
        "D": -3.5,
        "E": -3.5,
        "K": -3.9,
        "R": -4.5,
        "H": -3.2,
    }

    seq_len = len(AA_seq_segment)
    edge_weight = Edge_weight
    center_weight = 100

    if Model == "Linear Variation":
        # Linear weight variation
        weights = np.concatenate((np.linspace(edge_weight, center_weight, num=seq_len // 2 + 1),
                                  np.linspace(center_weight, edge_weight, num=seq_len // 2 + 1)[1:]),
                                 axis=None)

        hydropathy_index = 0
        for index, AA in enumerate(AA_seq_segment):
            hydropathy_index += (Kyte_Doolittle_scale[AA] * weights[index])

        avg_hydropathy_index = round(hydropathy_index / np.sum(weights), 3)
        return avg_hydropathy_index
    else:
        # Exponential weight variation
        weights = np.concatenate((np.geomspace(edge_weight, center_weight, num=seq_len // 2 + 1),
                                  np.geomspace(center_weight, edge_weight, num=seq_len // 2 + 1)[1:]),
                                 axis=None)

        hydropathy_index = 0
        for index, AA in enumerate(AA_seq_segment):
            hydropathy_index += (Kyte_Doolittle_scale[AA] * weights[index])

        avg_hydropathy_index = round(hydropathy_index / np.sum(weights), 3)
        return avg_hydropathy_index


def Hydropathicity_array_gen(AA_seq: str, Window_size: int, EDGE_weight=100, model="Linear Variation") -> (
        list[float], list[int]):
    hyrdopathicity_array = []
    AA_range = np.arange(start=(Window_size // 2),
                         stop=(len(AA_seq) - (Window_size) // 2))
    for index in AA_range:
        hyrdopathicity_array.append(Hydropathicity_value_calc(
            AA_seq_segment=AA_seq[index - (Window_size // 2): index + (Window_size // 2 + 1)],
            Edge_weight=EDGE_weight,
            Model=model))

    return np.array(hyrdopathicity_array), AA_range


######################################################################################## FRONTEND #######################################################################

st.set_page_config(
   page_title="HYDROPATHY INDEX PLOTTER ",
   page_icon="📈",
   layout="wide")

Title = st.container()
Intro = st.container()
Input_box = st.form("input")
Outbut_box = st.container()
Plot = st.container()

st.markdown("<style>.main {background-color: #FFFFFF;color:black;} </style>", unsafe_allow_html=True)

with Title:
    Title.markdown(
        "<h1 style='font-family : century gothic;text-align: center; color: Black;'><u>HYDROPATHY INDEX PLOTTER<u></h1>",
        unsafe_allow_html=True)

with Input_box:
    Input_box.markdown(
        "<h4 style='font-family : century gothic;text-align: center; color: Black;'><u>ENTER YOUR SEQUENCE or ACCESSION ID<u></h4>",
        unsafe_allow_html=True)
    input = Input_box.text_area(
        label="Enter NCBI Accession ID of your peptide sequence or Paste the FASTA sequence in the box",
        value="NP_001035835.1")
    input_type = Input_box.radio(label="| Select Your Input Type |", options=["ACCESSION ID", "FASTA"])
    computation_model = Input_box.selectbox(label="| Computation Model |",
                                            options=["Linear Variation", "Exponential Variation"])
    window_size = Input_box.slider(label="| Window Size |", min_value=3, max_value=21, step=2, value=7)
    edge_weight = Input_box.slider(label="| Edge Weight |", min_value=1, max_value=100, step=1, value=100)
    input_submit = st.form_submit_button("Submit")

    with st.expander("*Click here to know more about the parameters"):
        st.markdown("<h6 style='font-family : Garamond;text-align: center; color: #778899;'><u>To be Updated<u></h6>",unsafe_allow_html=True)

    if input_submit:
        if input_type == "FASTA":
            chk = FASTA_Parser(input.strip())
            if chk == -1:
                Input_box.error("Please, Enter a Valid FASTA Sequence")
            else:
                SEQ = chk.upper()
        else:
            chk = Fetch_Sequence_NCBI(input.strip())
            if chk == -1:
                Input_box.error("Please, Enter a Valid Acession ID")
            else:
                SEQ = FASTA_Parser(chk).upper()

        if chk != -1:
            Outbut_box.markdown(
                "<h4 style='font-family : century gothic;text-align: center; color: Black;'><u>PEPTIDE SEQUENCE<u></h4>",
                unsafe_allow_html=True)
            Outbut_box.text(Preview_Sequence(SEQ))

with Plot:
    if input_submit == True and chk != -1:
        Plot.markdown(
            "<h4 style='font-family : century gothic;text-align: center; color: Black;'><u>HYDROPATHY PLOT<u></h4>",
            unsafe_allow_html=True)

        try:
            hydropathicity_values, Amino_Acid_range = Hydropathicity_array_gen(SEQ,
                                                                               Window_size=window_size,
                                                                               EDGE_weight=edge_weight,
                                                                               model=computation_model)
            validity = True
        except KeyError:
            validity = False
            Plot.error("Invalid Amino Acid code present")

        if validity:
            plt.style.use("ggplot")
            fig, ax = plt.subplots()
            fig.set_size_inches(18.5, 12)
            fig.set_dpi(800)
            ax.hlines(y=0,
                      xmin=0, xmax=len(SEQ),
                      color="black",
                      linestyles="--",
                      alpha=0.5)
            ax.plot(Amino_Acid_range, hydropathicity_values,
                    color="black",
                    alpha=0.75)

            ax.set_xlim(0, len(SEQ))
            ax.xaxis.set_minor_locator(MultipleLocator(5))
            ax.yaxis.set_minor_locator(MultipleLocator(0.2))
            ax.tick_params(bottom=True,
                           top=True,
                           left=True,
                           right=True)
            ax.tick_params(labelbottom=True,
                           labeltop=True,
                           labelleft=True,
                           labelright=True)

            ax.fill_between(Amino_Acid_range, hydropathicity_values,
                            where=(hydropathicity_values >= 0),
                            interpolate=True,
                            alpha=0.75,
                            label="Hydrophobic Regions")
            ax.fill_between(Amino_Acid_range, hydropathicity_values,
                            where=(hydropathicity_values <= 0),
                            color="teal",
                            interpolate=True,
                            alpha=0.75,
                            label="Hydrophilic   Regions")

            ax.set_xlabel("Amino Acid Number →", fontsize="18")
            ax.set_ylabel("← Hydropathy Index →", fontsize="18")
            ax.legend(fontsize=15, edgecolor="black", facecolor="white")

            Plot.pyplot(fig)

            Plot.markdown("***")
            Plot.markdown("Check out my other projects at my [github](https://github.com/baby-phage).")
            Plot.markdown("***")
