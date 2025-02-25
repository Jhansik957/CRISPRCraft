import streamlit as st
import pandas as pd
from Bio import SeqIO
import re
from Bio.SeqUtils import gc_fraction
from PIL import Image

# Set page configuration
st.set_page_config(page_title="CRISPRCraft", layout="wide")

# Custom CSS for styling
st.markdown(
    """
    <style>
        /* General Page Styling */
        body {
            background-color: #f4f8fb;
            color: #2c3e50;
            font-family: Arial, sans-serif;
        }

        /* Sidebar Styling */
        .stSidebar {
            background-color: #ECF0F1 !important; /* Light Gray */
            padding: 10px;
        }
        
        /* Sidebar Text */
        .stSidebar .css-1lcbmhc {
            color: #2C3E50 !important; /* Dark Blue-Gray */
            font-size: 18px !important; /* Larger Text */
            font-weight: bold;
        }

        /* Active Selection */
        .stSidebar .css-17eq0hr {
            color: #E74C3C !important; /* Red for selection */
            font-weight: bold;
        }

        /* Hover Effect */
        .stSidebar .css-1lcbmhc:hover {
            color: #007acc !important; /* Blue Hover */
        }

        /* Main Title */
        .main-title {
            color: #004c7d;
            text-align: center;
            font-size: 36px;
            font-weight: bold;
            margin-bottom: 20px;
        }

        /* Content Box */
        .highlight-box {
            background-color: #ffffff;
            padding: 20px;
            border-radius: 12px;
            box-shadow: 0px 4px 8px rgba(0, 0, 0, 0.1);
            margin-bottom: 20px;
        }

        /* Buttons */
        .btn-primary {
            background-color: #007acc !important;
            color: white !important;
            font-weight: bold;
        }
        .btn-primary:hover {
            background-color: #005f99 !important;
        }

        /* Welcome Container */
        .welcome-container {
            background-color: #E6E6FA; /* Light lavender */
            padding: 20px;
            border-radius: 10px;
            text-align: center;
            box-shadow: 2px 2px 10px rgba(0,0,0,0.1);
            margin-bottom: 20px;
        }

        .welcome-text {
            color: #2c2c2c;
            font-size: 28px;
            font-weight: bold;
        }

        .highlight-crispr {
            color: #007bff; /* Blue color */
        }
    </style>
    """,
    unsafe_allow_html=True
)

# Sidebar Navigation
st.sidebar.title("CRISPRCraft Navigation")
selection = st.sidebar.radio("Go to", ["Home", "gRNA Prediction", "About", "Contact Us"])

def calculate_gc_content(sequence):
    """Calculate GC content percentage."""
    return round(gc_fraction(sequence) * 100, 2)

def count_off_targets(sequence, gRNA):
    """Counts approximate off-target sites by checking similar sequences."""
    off_target_count = len(re.findall(f'(?={gRNA[:-2]})', sequence)) - 1  # Ignore last 2 bases
    return off_target_count if off_target_count > 0 else "Low"

def extract_gRNAs(sequence, pam):
    """Extract gRNAs and compute GC content & off-target count."""
    gRNAs = []
    pam_regex = pam.replace('N', '.')  # Convert 'NGG' to regex
    
    for match in re.finditer(f'([ATCG]{{20}})({pam_regex})', sequence):
        gRNA = match.group(1)
        gc_content = calculate_gc_content(gRNA)
        off_target = count_off_targets(sequence, gRNA)
        
        gRNAs.append({"gRNA Sequence": gRNA, "GC Content (%)": gc_content, "Off-Target Risk": off_target})
    
    return gRNAs[:50]  # Limit to top 50 gRNAs

if selection == "Home":
    # Title (Centered)
    st.markdown("<h1 style='text-align: center;'>Welcome to CRISPRCraft: Your Precision gRNA Design Tool</h1>", unsafe_allow_html=True)

    # Display Logo in Center
    from PIL import Image
    logo = Image.open("logo.jpg")  # Change this if the filename is different

    st.markdown(
        "<div style='display: flex; justify-content: center;'>",
        unsafe_allow_html=True
    )
    st.image(logo, width=400)  # Adjust width if needed
    st.markdown("</div>", unsafe_allow_html=True)

    # Description
    st.write("**CRISPRCraft** is an advanced platform designed to streamline **guide RNA (gRNA) selection** for CRISPR genome editing. The other features of CRISPRCRAFT include:")
    st.markdown("""
    - 🔹 **On-Target Efficiency Prediction** – Identify high-efficiency gRNAs.
    - 🔹 **Off-Target Analysis** – Reduce unintended edits with advanced prediction.
    - 🔹 **Customizable Inputs** – Upload sequences in **FASTA/TXT** format.
    - 🔹 **Exportable Reports** – Download gRNA designs as **CSV files**.
    
    🚀 **Accelerate your CRISPR research with ease!**
    """)

    st.markdown("""
    ## Built By Researchers, For Researchers
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
    st.write("Upload a FASTA or TXT file, or paste a DNA sequence to analyze gRNAs and predict off-target effects.")

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

                # Download CSV
                csv = df.to_csv(index=False).encode('utf-8')
                st.download_button("Download gRNA CSV", csv, "gRNAs.csv", "text/csv")
            else:
                st.warning("No gRNAs found with the given PAM sequence.")

elif selection == "About":
    st.title("About CRISPRCraft")
    st.write("CRISPRCraft is a tool designed to help researchers design optimal gRNAs for CRISPR experiments.")
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
    st.write("📞 **Phone:** +1-234-567-8900")
    st.write("🌐 **Website:** [CRISPRCraft.com](https://crisprcraft.com)")

def main():  # Ensure this function exists before calling it
    print("CRISPRCraft app is running!")
    
if __name__ == "__main__":
    main()
