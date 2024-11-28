from flask import Flask, render_template, request, send_file
from Bio import Entrez
import pandas as pd
import os

app = Flask(__name__)
app.secret_key = 'supersecretkey'

# Set your email for PubMed API
Entrez.email = "your_email@example.com"

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/search', methods=['POST'])
def search():
    # Get the department query from the form
    query = request.form.get('query')

    # Fetch data from PubMed
    handle = Entrez.esearch(db="pubmed", term=query, retmax=50, sort="pub_date")
    record = Entrez.read(handle)
    handle.close()

    ids = record["IdList"]
    results = []

    for pub_id in ids:
        summary_handle = Entrez.esummary(db="pubmed", id=pub_id)
        summary = Entrez.read(summary_handle)
        summary_handle.close()

        for doc in summary:
            results.append({
                "DOI": doc.get("DOI", "N/A"),
                "PMID": doc["Id"],
                "Title": doc["Title"],
                "Journal": doc["Source"],
                "Year": doc["PubDate"].split(" ")[0] if doc.get("PubDate") else "N/A"
            })

    # Convert to DataFrame
    df = pd.DataFrame(results)

    # Save as CSV
    output_csv = "results.csv"
    df.to_csv(output_csv, index=False)

    # Save as text
    output_text = "results.txt"
    df.to_string(open(output_text, "w"), index=False)

    return render_template('results.html', results=results, csv=output_csv, text=output_text)

@app.route('/download/<filename>')
def download(filename):
    return send_file(filename, as_attachment=True)

if __name__ == '__main__':
    app.run(debug=True)
