import pandas as pd
import os
from config import *


def generate_html_table(df, output_file='output.html'):
    """
    Generate an interactive HTML table from a pandas DataFrame
    """
    # Print DataFrame for debugging
    print("DataFrame contents:")
    print(df)
    print("\nDataFrame columns:")
    print(df.columns)

    # Generate table headers first
    headers = ""
    for i, col in enumerate(df.columns):
        headers += f'<th onclick="sortTable({i})">{col}</th>\n'

    # Generate table rows
    rows = ""
    for index, row in df.iterrows():
        row_html = "<tr>"

        # First column (fold pair) with link
        fold_pair = str(row.iloc[0])
        row_html += f'<td><a href="https://steveabecassis.github.io/MsaCluster/HTML/{fold_pair}.html" target="_blank">{fold_pair}</a></td>'

        # Remaining columns
        for value in row.iloc[1:]:
            try:
                # Try to convert to float and check if it's a number
                float_val = float(value)
                if pd.notna(float_val):
                    cell_value = f"{float_val:.2f}"
                else:
                    cell_value = "-"
            except (ValueError, TypeError):
                # If conversion fails, check if it's a string or other non-null value
                if pd.notna(value):
                    cell_value = str(value)
                else:
                    cell_value = "-"
            row_html += f"<td>{cell_value}</td>"

        row_html += "</tr>\n"
        rows += row_html

    # Create the complete HTML content (rest of the HTML template remains the same)
    html_content = f"""
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <title>Interactive Protein Comparison Table</title>
        <style>
            body {{
                font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
                background-color: #121212;
                color: #E0E0E0;
                margin: 20px;
            }}
            table {{
                width: 80%;
                margin: auto;
                border-collapse: collapse;
                box-shadow: 0 4px 8px rgba(0,0,0,0.5);
            }}
            th, td {{
                border: 1px solid #333333;
                padding: 12px 16px;
                text-align: left;
            }}
            th {{
                background-color: #004d40;
                color: #ffffff;
                font-size: 16px;
                cursor: pointer;
            }}
            tr:nth-child(even) {{
                background-color: #424242;
            }}
            tr:hover {{
                background-color: #616161;
            }}
            a {{
                color: #BB86FC;
                text-decoration: none;
            }}
            a:hover {{
                text-decoration: underline;
            }}
            h2 {{
                text-align: center;
                color: #E0E0E0;
                margin-top: 0;
            }}
        </style>
    </head>
    <body>
        <h2>Results Table</h2>
        <table id="myTable">
            <thead>
                <tr>
                    {headers}
                </tr>
            </thead>
            <tbody>
                {rows}
            </tbody>
        </table>

        <script>
            var sortDirection = [];

            function sortTable(column) {{
                var table, rows, switching, i, x, y, xValue, yValue, shouldSwitch;
                table = document.getElementById("myTable");
                switching = true;
                sortDirection[column] = !sortDirection[column];

                while (switching) {{
                    switching = false;
                    rows = table.getElementsByTagName("TR");
                    for (i = 1; i < rows.length - 1; i++) {{
                        shouldSwitch = false;
                        x = rows[i].getElementsByTagName("TD")[column];
                        y = rows[i + 1].getElementsByTagName("TD")[column];
                        xValue = x.innerHTML === '-' ? -Infinity : (isNaN(parseFloat(x.innerHTML)) ? x.innerHTML.toLowerCase() : parseFloat(x.innerHTML));
                        yValue = y.innerHTML === '-' ? -Infinity : (isNaN(parseFloat(y.innerHTML)) ? y.innerHTML.toLowerCase() : parseFloat(y.innerHTML));

                        if (sortDirection[column]) {{
                            if (xValue > yValue) {{
                                shouldSwitch = true;
                                break;
                            }}
                        }} else {{
                            if (xValue < yValue) {{
                                shouldSwitch = true;
                                break;
                            }}
                        }}
                    }}
                    if (shouldSwitch) {{
                        rows[i].parentNode.insertBefore(rows[i + 1], rows[i]);
                        switching = true;
                    }}
                }}
            }}
        </script>
    </body>
    </html>
    """

    # Write to file
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write(html_content)


def generate_html_table_from_parquet(summary_table_file, similarity_table_file, output_file):
    df = pd.read_parquet(summary_table_file)
    # Add the column of similarity between the two fold of the fold pair (score calulated with the script get_tm_align_score.py)
    fold1_fold2_sim = pd.read_parquet(similarity_table_file)
    df_all = pd.merge(df, fold1_fold2_sim, on='fold_pair')
    # Generate the HTML file
    generate_html_table(df_all, output_file)


def main():
    # Dataframe summary results generated by the script summary_table.py
    # Set output file
    output_file = MAIN_DIR + "/protein_comparison_table.html"

    print(f"Reading parquet file and generating")
    generate_html_table_from_parquet(SUMMARY_RESULTS_TABLE, SIMILARITY_RESULTS_TABLE, output_file)
    print("HTML file " + output_file + " generated successfully!")


if __name__ == "__main__":
    main()


