from flask import Flask, render_template, request, send_file
import os
import pandas as pd
from Bio import Entrez

app = Flask(__name__)  # Define the Flask application instance

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/search', methods=['POST'])
def search():
    query = request.form['query']
    email = request.form['email']
    Entrez.email = email  # Update email dynamically based on user input

    # Search parameters
    results_per_request = 1000
    total_results_to_fetch = 10000

    # PubMed search logic...
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=total_results_to_fetch, usehistory="y")
        record = Entrez.read(handle)
        handle.close()

        id_list = record["IdList"]

        if not id_list:
            return render_template('results.html', message="No results found.", total_files_created=0)

        os.makedirs('downloads', exist_ok=True)
        combined_data = []
        total_files = 0
        individual_files = []

        for i in range(0, len(id_list), results_per_request):
            chunk_ids = id_list[i:i + results_per_request]

            summary_handle = Entrez.esummary(db="pubmed", id=",".join(chunk_ids), retmax=results_per_request)
            summaries = Entrez.read(summary_handle)
            summary_handle.close()

            data = [
                {
                    'PMID': doc.get('Id', ''),
                    'Title': doc.get('Title', ''),
                    'Source': doc.get('Source', ''),
                    'PubDate': doc.get('PubDate', '')
                } for doc in summaries
            ]

            df = pd.DataFrame(data)
            total_files += 1

            csv_filename = f'results_{total_files}.csv'
            text_filename = f'results_{total_files}.txt'
            df.to_csv(f'downloads/{csv_filename}', index=False)
            df.to_csv(f'downloads/{text_filename}', index=False, sep='\t')

            individual_files.append(csv_filename)
            individual_files.append(text_filename)

            combined_data.extend(data)

        combined_df = pd.DataFrame(combined_data)
        combined_csv_path = 'downloads/combined_results.csv'
        combined_txt_path = 'downloads/combined_results.txt'
        combined_df.to_csv(combined_csv_path, index=False)
        combined_df.to_csv(combined_txt_path, index=False, sep='\t')

        return render_template(
            'results.html',
            total_files_created=total_files,
            combined_csv="combined_results.csv",
            combined_txt="combined_results.txt",
            individual_files=individual_files
        )
    except Exception as e:
        return render_template('results.html', message=f"An error occurred: {str(e)}", total_files_created=0)

@app.route('/download/<filename>')
def download(filename):
    file_path = os.path.join('downloads', filename)
    if os.path.exists(file_path):
        return send_file(file_path, as_attachment=True)
    else:
        return "File not found", 404

if __name__ == '__main__':
    app.run(debug=True)
