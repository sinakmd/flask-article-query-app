from flask import Flask, render_template, request
import pandas as pd
import openpyxl
import os
from Bio import Entrez

app = Flask(__name__)

# Set up Entrez email
Entrez.email = "your-email@example.com"

@app.route("/")
def index():
    return render_template("index.html")

@app.route("/search", methods=["POST"])
def search():
    query = request.form.get("query")
    results = []

    if query:
        # Search PubMed using Entrez
        handle = Entrez.esearch(db="pubmed", term=query, retmax=10)
        record = Entrez.read(handle)
        handle.close()

        id_list = record.get("IdList", [])
        for pubmed_id in id_list:
            # Fetch details for each PubMed ID
            handle = Entrez.esummary(db="pubmed", id=pubmed_id)
            summary = Entrez.read(handle)
            handle.close()

            if summary and "DocumentSummarySet" in summary:
                doc = summary["DocumentSummarySet"]["DocumentSummary"][0]
                results.append({
                    "DOI": doc.get("ELocationID", ""),
                    "PMID": pubmed_id,
                    "Title": doc.get("Title", ""),
                    "Journal": doc.get("Source", ""),
                    "Year": doc.get("PubDate", "").split(" ")[0]
                })

    # Convert results to a DataFrame
    df = pd.DataFrame(results)

    # Save to CSV and Excel
    csv_file = "results.csv"
    excel_file = "results.xlsx"
    df.to_csv(csv_file, index=False)
    df.to_excel(excel_file, index=False, engine="openpyxl")

    return render_template("results.html", results=results)

if __name__ == "__main__":
    # Bind to 0.0.0.0 and use the PORT environment variable
    port = int(os.environ.get("PORT", 5000))  # Default to port 5000 if not set
    app.run(host="0.0.0.0", port=port)
