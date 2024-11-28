import os
from flask import Flask, render_template, request, jsonify
import requests
import pandas as pd
from io import StringIO
import codecs

app = Flask(__name__)

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/search', methods=['POST'])
def search():
    try:
        query = request.form.get('query')
        email = request.form.get('email')

        # Simulate fetching PubMed data
        results = [
            {"Title": "Example Title 1", "Authors": "Author1, Author2", "Journal": "Journal1", "Year": "2021"},
            {"Title": "Example Title 2", "Authors": "Author3, Author4", "Journal": "Journal2", "Year": "2022"},
            # Add more simulated data here
        ]

        # Save combined CSV and TXT
        output_dir = 'downloads'
        os.makedirs(output_dir, exist_ok=True)

        combined_csv = "combined_results.csv"
        combined_txt = "combined_results.txt"

        with codecs.open(os.path.join(output_dir, combined_csv), 'w', encoding='utf-8') as csv_file:
            pd.DataFrame(results).to_csv(csv_file, index=False)

        with codecs.open(os.path.join(output_dir, combined_txt), 'w', encoding='utf-8') as txt_file:
            pd.DataFrame(results).to_csv(txt_file, sep='\t', index=False)

        # Return results
        return render_template(
            'results.html',
            results=results,
            total_files_created=1,  # Simulated total files created
            combined_csv=combined_csv,
            combined_txt=combined_txt
        )

    except Exception as e:
        return f"An error occurred: {e}", 500

@app.route('/download/<filename>')
def download_file(filename):
    file_path = os.path.join('downloads', filename)
    if os.path.exists(file_path):
        return send_file(file_path, as_attachment=True)
    else:
        return f"File {filename} not found!", 404

if __name__ == '__main__':
    # Use the PORT environment variable set by the hosting platform
    port = int(os.environ.get('PORT', 5000))
    app.run(host='0.0.0.0', port=port, debug=True)
