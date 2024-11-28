from flask import Flask, render_template, request, send_file, jsonify
import os
import pandas as pd
from Bio import Entrez
import time

app = Flask(__name__)

# Configure NCBI Entrez
Entrez.email = None  # Placeholder email; replaced dynamically during runtime

# Ensure downloads directory exists
DOWNLOADS_DIR = 'downloads'
os.makedirs(DOWNLOADS_DIR, exist_ok=True)

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/search', methods=['POST'])
def search():
    query = request.form['query']
    email = request.form['email']
    Entrez.email = email  # Set user-provided email
    results_per_request = 1000  # Max results per request for PubMed API
    max_results = 10000  # PubMed's hard limit for total results

    # Fetch search results
    try:
        search_handle = Entrez.esearch(
            db="pubmed", term=query, retmax=max_results, usehistory="y"
        )
        search_results = Entrez.read(search_handle)
        search_handle.close()

        # Extract WebEnv and QueryKey for subsequent paginated requests
        webenv = search_results["WebEnv"]
        query_key = search_results["QueryKey"]
        count = int(search_results["Count"])

        # Calculate the total number of files needed
        total_files = (count + results_per_request - 1) // results_per_request
        total_files = min(total_files, max_results // results_per_request)

        # Fetch summaries in batches and save to files
        for i in range(total_files):
            start = i * results_per_request + 1
            handle = Entrez.esummary(
                db="pubmed", retstart=start - 1, retmax=results_per_request,
                webenv=webenv, query_key=query_key
            )
            records = Entrez.read(handle)
            handle.close()

            # Process records into a pandas DataFrame
            articles = []
            for record in records:
                articles.append({
                    "DOI": record.get("DOI", "N/A"),
                    "PMID": record["Id"],
                    "Title": record.get("Title", "N/A"),
                    "Journal": record.get("Source", "N/A"),
                    "Year": record.get("PubDate", "N/A").split()[0],
                })

            # Create CSV and text files
            df = pd.DataFrame(articles)
            csv_filename = f"{DOWNLOADS_DIR}/results_{i + 1}.csv"
            txt_filename = f"{DOWNLOADS_DIR}/results_{i + 1}.txt"
            df.to_csv(csv_filename, index=False)
            df.to_csv(txt_filename, index=False, sep='\t')

        return render_template(
            'results.html', total_files=total_files
        )

    except Exception as e:
        return f"An error occurred: {e}"

@app.route('/download/<file_type>/<filename>')
def download(file_type, filename):
    # Validate file type
    if file_type not in ['csv', 'txt']:
        return "Invalid file type!", 400

    # Construct the file path
    file_path = os.path.join(DOWNLOADS_DIR, filename)
    if os.path.exists(file_path):
        return send_file(file_path, as_attachment=True)
    else:
        return f"File {filename} not found!", 404

if __name__ == '__main__':
    app.run(debug=True)
