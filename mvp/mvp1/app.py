from flask import Flask, render_template, request, send_file
import os
import pandas as pd
import requests

app = Flask(__name__)

# Function to fetch PubMed results with pagination
def fetch_pubmed_results(query, email):
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
    retmax = 1000  # Maximum results per request
    max_results = 10000  # API hard limit
    results = []

    search_url = f"{base_url}esearch.fcgi"
    params = {
        "db": "pubmed",
        "term": query,
        "retmode": "json",
        "retmax": retmax,
        "email": email
    }

    retstart = 0  # Start from the first result

    # Fetch results in paginated chunks
    while retstart < max_results:
        params["retstart"] = retstart
        response = requests.get(search_url, params=params)
        response.raise_for_status()
        search_data = response.json()

        # Extract IDs
        id_list = search_data.get("esearchresult", {}).get("idlist", [])
        total_results = int(search_data.get("esearchresult", {}).get("count", 0))

        # Append IDs to results
        results.extend(id_list)

        # Break if we have fetched all available results
        if len(results) >= total_results or len(id_list) == 0:
            break

        retstart += retmax  # Move to the next batch

    # Limit total results to the API's hard limit
    results = results[:max_results]
    print(f"Total results fetched: {len(results)}")

    # Fetch summaries in chunks
    summaries = []
    chunk_size = 200
    for i in range(0, len(results), chunk_size):
        chunk = results[i:i + chunk_size]
        summary_url = f"{base_url}esummary.fcgi"
        summary_params = {
            "db": "pubmed",
            "id": ",".join(chunk),
            "retmode": "json",
            "email": email
        }
        try:
            summary_response = requests.get(summary_url, params=summary_params)
            summary_response.raise_for_status()
            chunk_summaries = summary_response.json().get("result", {})
            summaries.extend(chunk_summaries.values())
        except requests.exceptions.RequestException as e:
            print(f"Error during summary API call: {e}")
            continue

    # Format results
    formatted_results = [
        {
            "DOI": article.get("elocationid", ""),
            "PMID": article.get("uid", ""),
            "Title": article.get("title", ""),
            "Journal": article.get("source", ""),
            "Year": article.get("pubdate", "").split(" ")[0]
        }
        for article in summaries if isinstance(article, dict)
    ]

    return formatted_results


@app.route('/')
def index():
    return render_template('index.html')


@app.route('/search', methods=['POST'])
def search():
    query = request.form['query']
    email = request.form['email']

    print(f"Starting search for query: {query}")
    print(f"Email used: {email}")

    # Fetch PubMed results
    results = fetch_pubmed_results(query, email)

    if not results:
        return render_template('results.html', results=[], total_files=0)

    # Save results to CSV files
    os.makedirs('downloads', exist_ok=True)
    file_count = 0
    chunk_size = 1000  # Split results into 1,000 per file
    for i in range(0, len(results), chunk_size):
        chunk = results[i:i + chunk_size]
        df = pd.DataFrame(chunk)

        # Save each chunk to a new file
        csv_filename = f"downloads/results_{file_count + 1}.csv"
        txt_filename = f"downloads/results_{file_count + 1}.txt"
        df.to_csv(csv_filename, index=False)
        df.to_csv(txt_filename, index=False, sep='\t')
        file_count += 1

    return render_template('results.html', results=results, total_files=file_count)


@app.route('/download/<filename>')
def download_file(filename):
    file_path = os.path.join('downloads', filename)
    if os.path.exists(file_path):
        return send_file(file_path, as_attachment=True)
    else:
        return f"File {filename} not found!", 404


if __name__ == '__main__':
    app.run(debug=True)
