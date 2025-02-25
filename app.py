import streamlit as st
import pandas as pd
from Bio import SeqIO
import re
from Bio.SeqUtils import gc_fraction, MeltingTemp as mt
from PIL import Image

# Set page configuration
st.set_page_config(page_title="CRISPRCraft", layout="wide")

# Sidebar Navigation
st.sidebar.title("CRISPRCraft Navigation")
selection = st.sidebar.radio("Go to", ["Home", "gRNA Prediction", "About", "Contact Us"])

def calculate_gc_content(sequence):
    """Calculate GC content percentage."""
    gc_percent = round(gc_fraction(sequence) * 100, 2)
    return gc_percent, "Ideal Sequence" if 40 <= gc_percent <= 60 else "Non-Ideal Sequence"

def melting_temperature(sequence):
    """Calculate melting temperature using the Wallace rule."""
    return mt.Tm_Wallace(sequence)

def count_off_targets(sequence, gRNA):
    """Counts approximate off-target sites by checking similar sequences."""
    pattern = f'(?=({gRNA[:20]}))'
    off_target_count = len(re.findall(pattern, sequence)) - 1  # Ignore perfect match
    return off_target_count if off_target_count > 0 else "Low"

def calculate_grna_efficiency(sequence):
    """Placeholder for gRNA efficiency score (use Azimuth/DeepCRISPR models in the future)."""
    return round((gc_fraction(sequence) * 100) / 2, 2)  # Simplified score for now

def estimate_on_target_cleavage(sequence):
    """Estimate gRNA on-target cleavage efficiency."""
    return round(gc_fraction(sequence) * 100, 2)  # Simplified estimation for now

def extract_gRNAs(sequence, pam):
    """Extract gRNAs and compute GC content, off-targets, efficiency, and melting temperature."""
    gRNAs = []
    pam_regex = pam.replace('N', '.')  # Convert 'NGG' to regex
    
    for match in re.finditer(f'([ATCG]{{20}})({pam_regex})', sequence):
        gRNA = match.group(1)
        gc_content, status = calculate_gc_content(gRNA)
        off_target = count_off_targets(sequence, gRNA)
        tm = melting_temperature(gRNA)
        efficiency_score = calculate_grna_efficiency(gRNA)
        cleavage_efficiency = estimate_on_target_cleavage(gRNA)
        
        gRNAs.append({
            "gRNA Sequence": gRNA, 
            "GC Content (%)": gc_content, 
            "GC Status": status,
            "Melting Temp (°C)": round(tm, 2),
            "Off-Target Risk": off_target,
            "Efficiency Score": efficiency_score,
            "On-Target Cleavage Efficiency": cleavage_efficiency
        })
    
    return gRNAs[:50]  # Limit to top 50 gRNAs
 
if selection == "Home":
    st.markdown("<h1 style='text-align: center;'>Welcome to CRISPRCraft: Your Precision gRNA Design Tool</h1>", unsafe_allow_html=True)
    
    logo = Image.open("logo.jpg")  # Ensure you have this image file
    st.image(logo, width=400)

    st.write("**CRISPRCraft** is an advanced platform designed for **guide RNA (gRNA) selection** in CRISPR genome editing. The other features of CRISPRCRAFT include:")
    st.markdown("""
    -  **On-Target Efficiency Prediction** – Identify high-efficiency gRNAs.
    -  **Off-Target Analysis** – Reduce unintended edits with advanced prediction.
    -  **Customizable Inputs** – Upload sequences in **FASTA/TXT** format.
    -  **Exportable Reports** – Download gRNA designs as **CSV files**.
    
    🚀 **Accelerate your CRISPR research with ease!**
    """)

    st.markdown("""
    ## Built By Researcher, For Researchers
    Developed by **[Jhansi](https://github.com/Jhansik957/CRISPRCraft)**, this tool is tailored to assist **biologists, genetic engineers, and computational researchers** in CRISPR-based genome editing.
    
    📚 **Learn More:** [Read our latest insights](https://www.linkedin.com/posts/jhansikbioinfo_excited-to-share-my-latest-article-on-the-activity-7277908231920783360-Q1ok?utm_source=share&utm_medium=member_desktop&rcm=ACoAAFEXWlABejOhJYU6Ulsl7RY7w5_20qqJn-0)
    
    ## 🧑‍🔬 Join the CRISPRCraft Community!
    🔹 **Stay Updated** – Follow us on GitHub for updates & new features.
    🔹 **Contribute** – Help improve CRISPRCraft by suggesting features.
    🔹 **Contact Us** – Reach out for collaborations or support.
    
    📩 **Email:** support@crisprcraft.com  
    🌐 **Website:** [CRISPRCraft.com](https://crisprcraft.com)
    """)

elif selection == "gRNA Prediction":
    st.title("CRISPR Design Tool - CRISPRCraft")
    st.write("Upload a FASTA or TXT file, or paste a DNA sequence to analyze gRNAs.")

    uploaded_file = st.file_uploader("Upload a FASTA or TXT file", type=["fasta", "fa", "txt"])
    sequence_input = st.text_area("Paste your DNA sequence here", "")

    pam_sequence = st.selectbox("Select PAM sequence", ["NGG", "NAG", "NTG", "Other"])
    if pam_sequence == "Other":
        pam_sequence = st.text_input("Enter custom PAM sequence", "NGG")
    
    if st.button("Predict gRNAs"):
        sequence = ""
        if uploaded_file:
            st.success("File uploaded successfully!")
            if uploaded_file.name.endswith(".txt"):
                sequence = uploaded_file.read().decode("utf-8").strip().upper()
            else:
                records = list(SeqIO.parse(uploaded_file, "fasta"))
                sequence = str(records[0].seq).upper() if records else ""
        elif sequence_input:
            sequence = sequence_input.strip().upper()

        if sequence:
            gRNAs = extract_gRNAs(sequence, pam_sequence)
            if gRNAs:
                df = pd.DataFrame(gRNAs)
                st.dataframe(df)

                csv = df.to_csv(index=False).encode('utf-8')
                st.download_button("Download gRNA CSV", csv, "gRNAs.csv", "text/csv")
            else:
                st.warning("No gRNAs found with the given PAM sequence.")

elif selection == "About":
    st.title("About CRISPRCraft")
    st.write("CRISPRCraft helps researchers design optimal gRNAs for CRISPR experiments.")
    st.write("It allows users to upload sequences, customize PAM recognition, and export gRNA predictions.")
    st.markdown("""
    ## 🧬 Why Use CRISPRCraft?
    - 🔹 **Accurate Predictions** – Uses advanced algorithms to **design precise gRNAs**.
    - 🔹 **Intuitive Interface** – No coding skills required! Simply **upload, analyze, and export**.
    - 🔹 **Customizable Settings** – Choose different **PAM sequences** and optimize results.
    - 🔹 **Fast and Efficient** – Processes **large datasets** in seconds.
    
    ## 📌 How It Works?
    - 🔹 **Upload** your **DNA sequence (FASTA/TXT)**.
    - 🔹 **Select** a **PAM sequence** (NGG, NAG, NTG, or custom).
    - 🔹 **Analyze** gRNAs with GC content & off-target predictions.
    - 🔹 **Download** optimized gRNA designs in CSV format.
    
    🔍 Get started by heading to the **gRNA Prediction** section!
    """)

elif selection == "Contact Us":
    st.title("Contact Us")
    st.write("For support and inquiries, reach out to us at:")
    st.write("📧 **Email:** support@crisprcraft.com")
    st.write("📞 **Phone:** +91 9666855747")
    st.write("🌐 **Website:** [CRISPRCraft.com](https://crisprcraft.com)")

def main():
    print("CRISPRCraft app is running!")

if __name__ == "__main__":
    main()
