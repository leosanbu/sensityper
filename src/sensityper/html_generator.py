#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
HTML Generator Module for Sensiscript/Sensityper
Generates interactive HTML reports from TSV output files
"""

import csv
import os
from typing import List, Dict, Optional


def _read_tsv_to_dict(tsv_path: str) -> List[Dict[str, str]]:
    """
    Reads a TSV file and returns a list of dictionaries.

    Args:
        tsv_path: Path to the TSV file

    Returns:
        List of dictionaries, one per row (excluding header)
    """
    data = []
    with open(tsv_path, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            data.append(row)
    return data


def _check_alert_content(alert_tsv_path: Optional[str]) -> List[Dict[str, str]]:
    """
    Checks if alert_output.tsv has data beyond header.

    Args:
        alert_tsv_path: Path to alert_output.tsv file (optional)

    Returns:
        List of alert rows if any exist, empty list otherwise
    """
    if not alert_tsv_path or not os.path.exists(alert_tsv_path):
        return []

    try:
        data = _read_tsv_to_dict(alert_tsv_path)
        # Filter out rows where Alert column is empty or just whitespace
        alerts = [row for row in data if row.get('Alert', '').strip()]
        return alerts
    except Exception:
        return []


def _extract_antibiotics_from_tsv(tsv_path: str) -> List[str]:
    """
    Extracts antibiotic names from TSV column headers.
    Looks for columns ending with '_NWT' or '_WT' and extracts the antibiotic name.

    Args:
        tsv_path: Path to the TSV file

    Returns:
        List of antibiotic names found in the TSV (in order of appearance)
    """
    antibiotics = []
    with open(tsv_path, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f, delimiter='\t')
        headers = reader.fieldnames

        for header in headers:
            if header.endswith('_NWT'):
                # Extract antibiotic name (everything before '_NWT')
                abx = header[:-4]  # Remove '_NWT'
                if abx not in antibiotics:
                    antibiotics.append(abx)

    return antibiotics


def _generate_alert_banner_html(alert_data: List[Dict[str, str]]) -> str:
    """
    Generates HTML for alert banner displaying isolates with XDR or None flags.

    Args:
        alert_data: List of alert rows from alert_output.tsv

    Returns:
        HTML string for alert banner
    """
    if not alert_data:
        # Show banner with "No clinical alerts" message
        alerts_html = '            <li>No clinical alerts</li>'
    else:
        # Header message
        header = '            <li style="list-style: none; font-weight: 400; margin-bottom: 10px;">There are no effective treatment options among recommended regimens for the following isolates:</li>'

        # Build simplified alert list
        alert_items = [header]
        for row in alert_data:
            alert_type = row.get('Alert', '').strip()
            isolate = row.get('isolate', 'Unknown')

            if alert_type == 'XDR':
                message = f'Isolate <strong>{isolate}</strong>: XDR (Extensively Drug-Resistant)'
            elif alert_type == 'None':
                message = f'Isolate <strong>{isolate}</strong>: Unable to assign treatment'
            else:
                message = f'Isolate <strong>{isolate}</strong>: {alert_type}'

            alert_items.append(f'            <li>{message}</li>')

        alerts_html = '\n'.join(alert_items)

    return f'''    <div class="alert-banner">
        <h3 onclick="toggleAlerts()">
            <span>⚠️ Clinical Alerts</span>
            <span class="alert-toggle-icon" id="alert-toggle">▼</span>
        </h3>
        <ul id="alert-list">
{alerts_html}
        </ul>
    </div>
'''


def _get_css_template() -> str:
    """
    Returns the complete CSS template including base styles, tab navigation, and alert banner.

    Returns:
        CSS string
    """
    return '''    <style>
    * { box-sizing: border-box; }

    body {
      font-family: Arial, sans-serif;
      font-size: 14px;
      margin: 0;
      padding: 40px 20px;
      background: #f8f9fa;
      min-height: 100vh;
    }

    .container {
      max-width: 98%;
      margin: 0 auto;
      background: white;
      border-radius: 2px;
      box-shadow: 0 1px 3px rgba(0,0,0,0.08);
      padding: 40px;
      border-top: 3px solid #2c3e50;
    }

    h2 {
      font-size: 28px;
      font-weight: 300;
      text-align: left;
      margin: 0 0 10px 0;
      color: #2c3e50;
      letter-spacing: -0.5px;
    }

    .subtitle {
      font-size: 14px;
      color: #7f8c8d;
      margin-bottom: 40px;
      font-weight: 400;
    }

    .controls {
      display: flex;
      align-items: center;
      gap: 24px;
      margin-bottom: 32px;
      padding: 24px;
      background: #fafbfc;
      border: 1px solid #e1e4e8;
      border-radius: 2px;
      flex-wrap: wrap;
    }

    .control-group {
      display: flex;
      flex-direction: column;
      gap: 6px;
    }

    .controls label {
      font-weight: 500;
      color: #24292e;
      font-size: 13px;
      text-transform: uppercase;
      letter-spacing: 0.5px;
    }

    .controls select {
      padding: 10px 14px;
      border: 1px solid #d1d5da;
      border-radius: 2px;
      font-size: 14px;
      background: white;
      color: #24292e;
      cursor: pointer;
      transition: all 0.15s ease;
      min-width: 200px;
      font-family: inherit;
    }

    .controls select:hover {
      border-color: #2c3e50;
    }

    .controls select:focus {
      outline: none;
      border-color: #2c3e50;
      box-shadow: 0 0 0 3px rgba(44, 62, 80, 0.1);
    }

    .controls input[type="checkbox"] {
      margin-right: 8px;
      cursor: pointer;
      width: 16px;
      height: 16px;
    }

    .controls label:has(input[type="checkbox"]) {
      display: inline-block;
      cursor: pointer;
      font-size: 14px;
      text-transform: none;
      letter-spacing: normal;
      padding: 8px 0;
      white-space: nowrap;
    }

    .controls label:has(input[type="checkbox"]):hover {
      color: #0366d6;
    }

    /* Toggle Switch Styles */
    .toggle-switch-container {
      text-align: right;
      margin: 12px 0;
    }

    .toggle-switch {
      display: inline-flex;
      align-items: center;
      gap: 10px;
      cursor: pointer;
      font-size: 14px;
      color: #24292e;
      font-weight: 500;
    }

    .toggle-switch input[type="checkbox"] {
      display: none;
    }

    .toggle-slider {
      position: relative;
      width: 44px;
      height: 24px;
      background-color: #d1d5da;
      border-radius: 24px;
      transition: background-color 0.3s ease;
      box-shadow: inset 0 1px 3px rgba(0, 0, 0, 0.1);
    }

    .toggle-slider::before {
      content: '';
      position: absolute;
      width: 18px;
      height: 18px;
      border-radius: 50%;
      background-color: white;
      top: 3px;
      left: 3px;
      transition: transform 0.3s ease;
      box-shadow: 0 2px 4px rgba(0, 0, 0, 0.2);
    }

    .toggle-switch input[type="checkbox"]:checked + .toggle-slider {
      background-color: #22863a;
    }

    .toggle-switch input[type="checkbox"]:checked + .toggle-slider::before {
      transform: translateX(20px);
    }

    .toggle-switch:hover .toggle-slider {
      background-color: #b8bdc3;
    }

    .toggle-switch input[type="checkbox"]:checked + .toggle-slider:hover {
      background-color: #1a6d2e;
    }

    .controls input {
      padding: 10px 14px;
      border: 1px solid #d1d5da;
      border-radius: 2px;
      font-size: 14px;
      flex: 1;
      min-width: 250px;
      transition: all 0.15s ease;
      font-family: inherit;
    }

    .controls input:focus {
      outline: none;
      border-color: #2c3e50;
      box-shadow: 0 0 0 3px rgba(44, 62, 80, 0.1);
    }

    .controls input::placeholder {
      color: #959da5;
    }

    .table-wrapper {
      overflow-x: auto;
      border: 1px solid #e1e4e8;
      border-radius: 2px;
    }

    table {
      border-collapse: collapse;
      width: 100%;
      background: white;
    }

    th, td {
      border: 1px solid #e1e4e8;
      padding: 12px 16px;
      vertical-align: top;
      font-size: 13px;
      text-align: left;
    }

    th {
      background: #fafbfc;
      color: #24292e;
      cursor: pointer;
      position: sticky;
      top: 0;
      font-size: 12px;
      font-weight: 600;
      text-transform: uppercase;
      letter-spacing: 0.5px;
      z-index: 10;
      border-bottom: 2px solid #d1d5da;
    }

    th:hover {
      background: #f3f4f6;
    }

    tr:nth-child(even) {
      background-color: #fafbfc;
    }

    tr:hover {
      background-color: #f6f8fa;
    }

    td[colname="treatment recommendation"] div {
      display: flex !important;
      flex-wrap: wrap !important;
      gap: 4px !important;
    }

    td[colname="treatment recommendation"] div span {
      padding: 4px 10px !important;
      font-size: 11px !important;
      font-weight: 500;
      white-space: nowrap;
      border-radius: 2px;
      flex-shrink: 0;
      text-transform: capitalize;
      letter-spacing: 0.3px;
    }

    td div.list-items {
      margin: 0 !important;
      padding-left: 0 !important;
      font-size: 13px !important;
      line-height: 1.4 !important;
    }

    td div.list-items div {
      margin: 0 !important;
      padding-left: 0 !important;
      line-height: 1.3 !important;
      white-space: nowrap;
      list-style: none !important;
      font-family: Arial, monospace;
      font-size: 13px;
    }

    td[colname$="R"]:not(:empty) {
      background-color: #fff5f5 !important;
    }

    td[colname$="S"]:not(:empty) {
      background-color: #f6f8fa !important;
    }

    td[colname="isolate"] {
      font-weight: 600;
      color: #24292e;
      font-family: Arial, monospace;
      font-size: 12px;
    }

    th[data-sort-dir="asc"]::after {
      content: " ▲";
      font-size: 10px;
    }

    th[data-sort-dir="desc"]::after {
      content: " ▼";
      font-size: 10px;
    }

    .table-wrapper::-webkit-scrollbar {
      height: 8px;
    }

    .table-wrapper::-webkit-scrollbar-track {
      background: #f6f8fa;
    }

    .table-wrapper::-webkit-scrollbar-thumb {
      background: #d1d5da;
      border-radius: 4px;
    }

    .table-wrapper::-webkit-scrollbar-thumb:hover {
      background: #959da5;
    }

    .info-badge {
      display: inline-block;
      padding: 2px 8px;
      background: #f6f8fa;
      border: 1px solid #e1e4e8;
      border-radius: 2px;
      font-size: 11px;
      color: #586069;
      margin-left: 8px;
      font-weight: 500;
    }

    /* Tab navigation styles */
    .tab-container {
      border-bottom: 2px solid #e1e4e8;
      margin-bottom: 20px;
    }

    .tab-nav {
      display: flex;
      gap: 10px;
      list-style: none;
      padding: 0;
      margin: 0;
    }

    .tab-nav button {
      background: none;
      border: none;
      padding: 12px 24px;
      cursor: pointer;
      font-size: 14px;
      font-weight: 500;
      color: #586069;
      border-bottom: 2px solid transparent;
      transition: all 0.2s;
    }

    .tab-nav button:hover {
      color: #24292e;
      border-bottom-color: #d1d5da;
    }

    .tab-nav button.active {
      color: #24292e;
      border-bottom-color: #0366d6;
    }

    .tab-content {
      display: none;
    }

    .tab-content.active {
      display: block;
    }

    /* Alert banner styles */
    .alert-banner {
      background-color: #fff3cd;
      border: 1px solid #ffc107;
      border-radius: 4px;
      padding: 15px 20px;
      margin-bottom: 20px;
      color: #856404;
    }

    .alert-banner h3 {
      margin-top: 0;
      color: #856404;
      font-size: 16px;
      cursor: pointer;
      user-select: none;
      display: flex;
      align-items: center;
      justify-content: space-between;
    }

    .alert-banner h3:hover {
      color: #664d03;
    }

    .alert-toggle-icon {
      font-size: 14px;
      transition: transform 0.2s ease;
      display: inline-block;
    }

    .alert-toggle-icon.collapsed {
      transform: rotate(-90deg);
    }

    .alert-banner ul {
      margin: 10px 0 0 0;
      padding-left: 20px;
      max-height: 400px;
      overflow-y: auto;
      overflow-x: hidden;
      transition: max-height 0.3s ease, opacity 0.3s ease, margin 0.3s ease;
      opacity: 1;
    }

    .alert-banner ul.collapsed {
      max-height: 0;
      opacity: 0;
      margin: 0;
      overflow: hidden;
    }

    .alert-banner ul::-webkit-scrollbar {
      width: 8px;
    }

    .alert-banner ul::-webkit-scrollbar-track {
      background: #ffeaa7;
      border-radius: 4px;
    }

    .alert-banner ul::-webkit-scrollbar-thumb {
      background: #fdcb6e;
      border-radius: 4px;
    }

    .alert-banner ul::-webkit-scrollbar-thumb:hover {
      background: #f39c12;
    }

    .alert-banner li {
      margin: 5px 0;
    }

    /* Treatment table specific styles */
    td[colname="recommended_1"],
    td[colname="recommended_2"] {
      font-weight: 600;
      color: #24292e;
    }

    td[colname="Recommended Treatment"],
    td[colname="Comment"],
    td[colname="Predicted Profile"] {
      white-space: normal;
      word-wrap: break-word;
      max-width: 300px;
    }

    td[colname="Predicted Profile"] div {
      display: flex !important;
      flex-wrap: wrap !important;
      gap: 4px !important;
    }

    td[colname="Predicted Profile"] div span {
      padding: 4px 10px !important;
      font-size: 11px !important;
      font-weight: 500;
      white-space: nowrap;
      border-radius: 2px;
      flex-shrink: 0;
      text-transform: capitalize;
      letter-spacing: 0.3px;
    }

    td.regimen-yes {
      background-color: #d4edda !important;
      font-weight: 600;
      color: #155724;
    }

    td.regimen-no {
      background-color: #f8f9fa !important;
      color: #6c757d;
    }

    /* Pagination styles */
    .pagination-controls {
      display: flex;
      align-items: center;
      gap: 16px;
      margin-top: 20px;
      padding: 16px;
      background: #fafbfc;
      border: 1px solid #e1e4e8;
      border-radius: 2px;
      justify-content: space-between;
    }

    .pagination-info {
      font-size: 13px;
      color: #586069;
    }

    .pagination-buttons {
      display: flex;
      gap: 8px;
    }

    .pagination-buttons button {
      padding: 6px 12px;
      border: 1px solid #d1d5da;
      border-radius: 2px;
      background: white;
      color: #24292e;
      cursor: pointer;
      font-size: 13px;
      transition: all 0.15s ease;
    }

    .pagination-buttons button:hover:not(:disabled) {
      background: #f6f8fa;
      border-color: #2c3e50;
    }

    .pagination-buttons button:disabled {
      opacity: 0.5;
      cursor: not-allowed;
    }

    .pagination-buttons button.active {
      background: #0366d6;
      color: white;
      border-color: #0366d6;
    }

    .page-size-selector {
      display: flex;
      align-items: center;
      gap: 8px;
    }

    .page-size-selector label {
      font-size: 13px;
      color: #24292e;
      font-weight: 500;
    }

    .page-size-selector select {
      padding: 6px 10px;
      border: 1px solid #d1d5da;
      border-radius: 2px;
      font-size: 13px;
      background: white;
      cursor: pointer;
    }
  </style>'''


def _get_javascript_template(table_id: str = 'resultsTable', antibiotics: Optional[List[str]] = None) -> str:
    """
    Returns the complete JavaScript template for table interactivity.

    Args:
        table_id: ID of the table element
        antibiotics: List of antibiotics for filtering (if applicable)

    Returns:
        JavaScript string
    """
    # Convert antibiotics list to JavaScript array string
    import json
    if antibiotics:
        abx_list_js = json.dumps(antibiotics)
    else:
        abx_list_js = "['ceftriaxone', 'azithromycin', 'ciprofloxacin', 'tetracycline', 'penicillin', 'spectinomycin', 'zoliflodacin', 'gepotidacin']"

    # Build antibiotic filter mapping if antibiotics provided
    abx_cols_js = '{}'
    if antibiotics:
        abx_mapping = {}
        for abx in antibiotics:
            abx_mapping[abx] = [f'{abx}R', f'{abx}S']

        # Convert to JavaScript object string
        abx_items = []
        for abx, cols in abx_mapping.items():
            cols_str = "['" + "', '".join(cols) + "']"
            abx_items.append(f"        '{abx}': {cols_str}")
        abx_cols_js = '{\n' + ',\n'.join(abx_items) + '\n      }'

    return f'''  <script>
    // Store original text content before transformation (separate for each table)
    let originalData = {{
      'resultsTable': [],
      'treatmentTable': []
    }};
    let currentPage = 1;
    let pageSize = 10;
    let currentTableId = '{table_id}';

    function sortTable(n) {{
      const table = document.getElementById(currentTableId);
      const tbody = table.querySelector('tbody');
      const rows = Array.from(tbody.querySelectorAll('tr'));

      // Determine current sort direction for this column
      const header = table.querySelectorAll('th')[n];
      const currentDir = header.getAttribute('data-sort-dir') || 'asc';
      const newDir = currentDir === 'asc' ? 'desc' : 'asc';

      // Clear all headers
      table.querySelectorAll('th').forEach(th => {{
        th.removeAttribute('data-sort-dir');
      }});

      // Set new direction on clicked header
      header.setAttribute('data-sort-dir', newDir);

      // Sort rows
      rows.sort((rowA, rowB) => {{
        const idxA = Array.from(tbody.querySelectorAll('tr')).indexOf(rowA);
        const idxB = Array.from(tbody.querySelectorAll('tr')).indexOf(rowB);

        const tableData = originalData[currentTableId] || [];
        const textA = (tableData[idxA] && tableData[idxA][n]) ? tableData[idxA][n].toLowerCase() : '';
        const textB = (tableData[idxB] && tableData[idxB][n]) ? tableData[idxB][n].toLowerCase() : '';

        if (newDir === 'asc') {{
          return textA.localeCompare(textB);
        }} else {{
          return textB.localeCompare(textA);
        }}
      }});

      // Update originalData array to match new order
      const tableData = originalData[currentTableId] || [];
      const newOriginalData = [];
      rows.forEach(row => {{
        const oldIdx = parseInt(row.getAttribute('data-original-idx'));
        newOriginalData.push(tableData[oldIdx]);
      }});

      // Re-append rows in sorted order and update indices
      rows.forEach((row, newIdx) => {{
        tbody.appendChild(row);
        row.setAttribute('data-original-idx', newIdx);
      }});

      // Update originalData
      originalData[currentTableId] = newOriginalData;

      // Reset to page 1 and update pagination after sorting
      currentPage = 1;
      updatePagination();
    }}

    function filterColumnsByAbx() {{
      const selected = Array.from(document.getElementById('show-cols-abx').selectedOptions)
                           .map(opt => opt.value.toLowerCase());
      const table = document.getElementById('resultsTable');
      const rows = table.getElementsByTagName('tr');
      const abxCols = {abx_cols_js};

      const headerRow = rows[0];
      const colNames = [];
      for (let c = 0; c < headerRow.cells.length; c++) {{
        colNames[c] = headerRow.cells[c].getAttribute('colname') || headerRow.cells[c].textContent;
      }}

      const visibleSet = new Set(['isolate', 'treatment recommendation']);

      if (selected.includes('all') || selected.length === 0) {{
        for (let i = 0; i < colNames.length; i++) {{
          visibleSet.add(colNames[i]);
        }}
      }} else {{
        selected.forEach(abx => {{
          if (abxCols[abx]) {{
            abxCols[abx].forEach(col => visibleSet.add(col));
          }}
        }});
      }}

      // Check WT toggle setting
      const showWT = document.getElementById('show-wt-cols').checked;

      for (let r = 0; r < rows.length; r++) {{
        const cells = rows[r].cells;
        for (let c = 0; c < cells.length; c++) {{
          const cname = colNames[c];
          const isWTColumn = cname && cname.endsWith('S') && cname !== 'isolate';

          if (visibleSet.has(cname)) {{
            // Show column if in visible set, but hide WT columns if toggle is off
            if (isWTColumn && !showWT) {{
              cells[c].style.display = 'none';
            }} else {{
              cells[c].style.display = '';
            }}
          }} else {{
            cells[c].style.display = 'none';
          }}
        }}
      }}
    }}

    function toggleWTColumns() {{
      // Re-apply antibiotic filter, which now respects the WT toggle setting
      filterColumnsByAbx();
    }}

    function searchTable() {{
      const input = document.getElementById('search-box-' + currentTableId);
      const filter = input.value.toLowerCase();
      const table = document.getElementById(currentTableId);
      const tbody = table.querySelector('tbody');
      const rows = tbody.querySelectorAll('tr');

      const tableData = originalData[currentTableId] || [];

      rows.forEach((row, idx) => {{
        let show = false;

        // Search in original data
        if (tableData[idx]) {{
          const rowText = tableData[idx].join(' ').toLowerCase();
          if (rowText.indexOf(filter) > -1) {{
            show = true;
          }}
        }}

        row.style.display = show ? '' : 'none';
      }});

      // Update pagination after search filtering
      updatePagination();
    }}

    function parsePredictedProfile(profileText) {{
      // Parse "ceftriaxone=YES, azithromycin=NO, ..." format
      const predictions = {{}};
      const parts = profileText.split(',').map(s => s.trim());

      parts.forEach(part => {{
        const match = part.match(/(.+?)=(.+)/);
        if (match) {{
          const abx = match[1].trim().toLowerCase();
          const value = match[2].trim().toUpperCase();
          predictions[abx] = value;
        }}
      }});

      return predictions;
    }}

    function createPredictedProfileBadges() {{
      const treatmentTable = document.getElementById('treatmentTable');
      if (!treatmentTable) return;

      const tbody = treatmentTable.querySelector('tbody');
      const rows = tbody.querySelectorAll('tr');

      // Store original data BEFORE transformation
      rows.forEach((row, rowIdx) => {{
        const cells = row.querySelectorAll('td');
        originalData['treatmentTable'][rowIdx] = [];
        row.setAttribute('data-original-idx', rowIdx);
        cells.forEach((cell, cellIdx) => {{
          originalData['treatmentTable'][rowIdx][cellIdx] = cell.textContent.trim();
        }});
      }});

      // Transform cells
      rows.forEach(row => {{
        const cells = row.querySelectorAll('td');

        // Find the Predicted Profile column (index 1 after isolate)
        cells.forEach((cell, idx) => {{
          if (cell.getAttribute('colname') === 'Predicted Profile') {{
            const profileText = cell.textContent.trim();
            if (!profileText) return;

            const predictions = parsePredictedProfile(profileText);
            const badgeContainer = document.createElement('div');
            badgeContainer.style.display = 'flex';
            badgeContainer.style.flexWrap = 'wrap';
            badgeContainer.style.gap = '4px';

            // Create badges for each antibiotic
            Object.keys(predictions).forEach(abx => {{
              const badge = document.createElement('span');
              badge.textContent = abx;
              badge.style.padding = '4px 10px';
              badge.style.borderRadius = '2px';
              badge.style.fontSize = '11px';
              badge.style.fontWeight = '500';
              badge.style.whiteSpace = 'nowrap';
              badge.style.letterSpacing = '0.3px';
              badge.style.textTransform = 'capitalize';

              if (predictions[abx] === 'YES') {{
                badge.style.backgroundColor = '#22863a';
                badge.style.color = 'white';
              }} else {{
                badge.style.backgroundColor = '#bcbcbc';
                badge.style.color = 'white';
              }}

              badgeContainer.appendChild(badge);
            }});

            cell.innerHTML = '';
            cell.appendChild(badgeContainer);
          }}
        }});
      }});
    }}

    function createRecommendationSubcells() {{
      const table = document.getElementById('resultsTable');
      if (!table) return;
      const tbody = table.querySelector('tbody');
      const rows = tbody.querySelectorAll('tr');
      const abxList = {abx_list_js};

      // Store original data BEFORE transformation
      rows.forEach((row, rowIdx) => {{
        const cells = row.querySelectorAll('td');
        originalData['resultsTable'][rowIdx] = [];
        row.setAttribute('data-original-idx', rowIdx);
        cells.forEach((cell, cellIdx) => {{
          originalData['resultsTable'][rowIdx][cellIdx] = cell.textContent.trim();
        }});
      }});

      // Transform treatment recommendation cells and mechanism cells
      rows.forEach((row, i) => {{
        const cells = row.querySelectorAll('td');
        if (cells.length < 2) return;

        const recCell = cells[1];
        const recText = recCell.textContent.toLowerCase().trim().replace(/,/g, ',');
        const recSet = new Set(recText.split(',').map(s => s.trim().replace(/^\\([^)]*\\)\\s*/, '')));

        const recContainer = document.createElement('div');
        recContainer.style.display = 'flex';
        recContainer.style.flexWrap = 'wrap';
        recContainer.style.gap = '4px';

        abxList.forEach(abx => {{
          const subcell = document.createElement('span');
          subcell.style.padding = '4px 10px';
          subcell.style.borderRadius = '2px';
          subcell.style.fontSize = '11px';
          subcell.style.fontWeight = '500';
          subcell.style.whiteSpace = 'nowrap';
          subcell.style.letterSpacing = '0.3px';
          subcell.textContent = abx;

          if (recSet.has(abx)) {{
            subcell.style.backgroundColor = '#22863a';
            subcell.style.color = 'white';
          }} else {{
            subcell.style.backgroundColor = '#bcbcbc';
            subcell.style.color = 'white';
          }}
          recContainer.appendChild(subcell);
        }});

        recCell.innerHTML = '';
        recCell.appendChild(recContainer);

        for (let c = 2; c < cells.length; c++) {{
          const td = cells[c];
          let content = td.textContent.trim();
          if (!content) continue;

          const items = content.split(/[\\/\\s]+/).filter(item => item.length > 0);
          if (items.length > 1) {{
            const listDiv = document.createElement('div');
            listDiv.className = 'list-items';

            items.forEach(item => {{
              const itemDiv = document.createElement('div');
              itemDiv.textContent = item;
              listDiv.appendChild(itemDiv);
            }});

            td.innerHTML = '';
            td.appendChild(listDiv);
          }}
        }}
      }});
    }}

    function switchTab(tabName) {{
      // Hide all tabs
      document.querySelectorAll('.tab-content').forEach(tab => {{
        tab.classList.remove('active');
      }});

      // Deactivate all tab buttons
      document.querySelectorAll('.tab-nav button').forEach(btn => {{
        btn.classList.remove('active');
      }});

      // Show selected tab
      document.getElementById(tabName + '-tab').classList.add('active');

      // Activate selected button
      event.target.classList.add('active');

      // Update current table ID for pagination
      if (tabName === 'resistance') {{
        currentTableId = 'resultsTable';
      }} else if (tabName === 'treatment') {{
        currentTableId = 'treatmentTable';
      }}

      // Reset pagination for the new table
      currentPage = 1;
      updatePagination();
    }}

    function updatePagination() {{
      const table = document.getElementById(currentTableId);
      if (!table) return;

      const tbody = table.querySelector('tbody');
      const rows = Array.from(tbody.querySelectorAll('tr'));

      // Get current search filter
      const searchInput = document.getElementById('search-box-' + currentTableId);
      const searchFilter = searchInput ? searchInput.value.toLowerCase() : '';

      // Build visibleRows by checking search matches against originalData
      // (independent of current display state to avoid stale pagination)
      const tableData = originalData[currentTableId] || [];
      const visibleRows = [];
      rows.forEach((row, idx) => {{
        let matchesSearch = true;

        if (searchFilter) {{
          matchesSearch = false;
          if (tableData[idx]) {{
            const rowText = tableData[idx].join(' ').toLowerCase();
            if (rowText.indexOf(searchFilter) > -1) {{
              matchesSearch = true;
            }}
          }}
        }}

        if (matchesSearch) {{
          visibleRows.push(row);
        }} else {{
          row.style.display = 'none';
        }}
      }});

      const totalRows = visibleRows.length;
      const totalPages = pageSize === -1 ? 1 : Math.ceil(totalRows / pageSize);

      // Adjust current page if needed
      if (currentPage > totalPages && totalPages > 0) {{
        currentPage = totalPages;
      }}
      if (currentPage < 1) {{
        currentPage = 1;
      }}

      // Calculate pagination indices
      const startIdx = (currentPage - 1) * pageSize;
      const endIdx = startIdx + pageSize;

      // Show/hide rows based on pagination
      if (pageSize === -1) {{
        // Show all rows
        visibleRows.forEach(row => {{
          row.style.display = '';
        }});
      }} else {{
        visibleRows.forEach((row, idx) => {{
          if (idx >= startIdx && idx < endIdx) {{
            row.style.display = '';
          }} else {{
            row.style.display = 'none';
          }}
        }});
      }}

      // Update pagination info
      const paginationInfo = document.getElementById('pagination-info-' + currentTableId);
      if (paginationInfo) {{
        if (pageSize === -1) {{
          paginationInfo.textContent = `Showing all ${{totalRows}} entries`;
        }} else {{
          const startEntry = totalRows === 0 ? 0 : startIdx + 1;
          const endEntry = Math.min(endIdx, totalRows);
          paginationInfo.textContent = `Showing ${{startEntry}} to ${{endEntry}} of ${{totalRows}} entries`;
        }}
      }}

      // Update page buttons (Previous/Next enable/disable)
      const prevBtn = document.getElementById('prev-page-' + currentTableId);
      const nextBtn = document.getElementById('next-page-' + currentTableId);

      if (prevBtn) prevBtn.disabled = currentPage === 1 || pageSize === -1;
      if (nextBtn) nextBtn.disabled = currentPage === totalPages || pageSize === -1;

      // Update page number buttons
      const paginationButtonsDiv = document.getElementById('pagination-buttons-' + currentTableId);
      const prevButton = document.getElementById('prev-page-' + currentTableId);
      const nextButton = document.getElementById('next-page-' + currentTableId);

      if (paginationButtonsDiv && pageSize !== -1) {{
        // Remove existing page number buttons (keep prev/next)
        const existingPageButtons = paginationButtonsDiv.querySelectorAll('button:not(#prev-page-' + currentTableId + '):not(#next-page-' + currentTableId + ')');
        existingPageButtons.forEach(btn => btn.remove());

        const maxButtons = 5;
        let startPage = Math.max(1, currentPage - Math.floor(maxButtons / 2));
        let endPage = Math.min(totalPages, startPage + maxButtons - 1);

        if (endPage - startPage < maxButtons - 1) {{
          startPage = Math.max(1, endPage - maxButtons + 1);
        }}

        // Insert page number buttons between prev and next
        for (let i = startPage; i <= endPage; i++) {{
          const btn = document.createElement('button');
          btn.textContent = i;
          btn.onclick = () => goToPage(i);
          if (i === currentPage) btn.classList.add('active');
          paginationButtonsDiv.insertBefore(btn, nextButton);
        }}
      }} else if (paginationButtonsDiv) {{
        // Remove page number buttons when showing all entries
        const existingPageButtons = paginationButtonsDiv.querySelectorAll('button:not(#prev-page-' + currentTableId + '):not(#next-page-' + currentTableId + ')');
        existingPageButtons.forEach(btn => btn.remove());
      }}
    }}

    function changePageSize(size) {{
      pageSize = parseInt(size);
      currentPage = 1;
      updatePagination();
    }}

    function goToPage(page) {{
      currentPage = page;
      updatePagination();
    }}

    function previousPage() {{
      const table = document.getElementById(currentTableId);
      if (!table) return;

      if (currentPage > 1) {{
        currentPage--;
        updatePagination();
      }}
    }}

    function nextPage() {{
      const table = document.getElementById(currentTableId);
      if (!table) return;

      const tbody = table.querySelector('tbody');
      const rows = Array.from(tbody.querySelectorAll('tr'));

      // Calculate visible rows based on search filter (not current display state)
      const searchInput = document.getElementById('search-box-' + currentTableId);
      const searchFilter = searchInput ? searchInput.value.toLowerCase() : '';

      const tableData = originalData[currentTableId] || [];
      let visibleCount = 0;
      rows.forEach((row, idx) => {{
        if (searchFilter) {{
          if (tableData[idx]) {{
            const rowText = tableData[idx].join(' ').toLowerCase();
            if (rowText.indexOf(searchFilter) > -1) {{
              visibleCount++;
            }}
          }}
        }} else {{
          visibleCount++;
        }}
      }});

      const totalPages = pageSize === -1 ? 1 : Math.ceil(visibleCount / pageSize);

      if (currentPage < totalPages) {{
        currentPage++;
        updatePagination();
      }}
    }}

    function toggleAlerts() {{
      const alertList = document.getElementById('alert-list');
      const toggleIcon = document.getElementById('alert-toggle');

      if (!alertList || !toggleIcon) return;

      if (alertList.classList.contains('collapsed')) {{
        alertList.classList.remove('collapsed');
        toggleIcon.classList.remove('collapsed');
        toggleIcon.textContent = '▼';
        localStorage.setItem('alertsCollapsed', 'false');
      }} else {{
        alertList.classList.add('collapsed');
        toggleIcon.classList.add('collapsed');
        toggleIcon.textContent = '▶';
        localStorage.setItem('alertsCollapsed', 'true');
      }}
    }}

    window.addEventListener('load', function() {{
      createRecommendationSubcells();
      createPredictedProfileBadges();
      updatePagination();

      // Restore alert collapse state from localStorage
      const isCollapsed = localStorage.getItem('alertsCollapsed') === 'true';
      if (isCollapsed) {{
        const alertList = document.getElementById('alert-list');
        const toggleIcon = document.getElementById('alert-toggle');
        if (alertList && toggleIcon) {{
          alertList.classList.add('collapsed');
          toggleIcon.classList.add('collapsed');
          toggleIcon.textContent = '▶';
        }}
      }}
    }});
  </script>'''


def _generate_resistance_table_html(data: List[Dict[str, str]], antibiotics: List[str]) -> str:
    """
    Generates the resistance profile table HTML.

    Args:
        data: List of dictionaries from TSV
        antibiotics: List of antibiotic names

    Returns:
        HTML string for resistance table
    """
    if not data:
        return '<p>No data available</p>'

    # Build table headers
    headers_html = []
    headers_html.append('            <th colname="isolate" onclick="sortTable(0)">Isolate</th>')
    # Same data as the treatment tab's "Predicted Profile" column, so use the same visible label.
    # colname stays "treatment recommendation" — it is the internal key used by the CSS and JS.
    headers_html.append('            <th colname="treatment recommendation" onclick="sortTable(1)">Predicted Profile</th>')

    col_idx = 2
    for abx in antibiotics:
        headers_html.append(f'            <th colname="{abx}R" onclick="sortTable({col_idx})">{abx.capitalize()} NWT</th>')
        col_idx += 1
        headers_html.append(f'            <th colname="{abx}S" onclick="sortTable({col_idx})">{abx.capitalize()} WT</th>')
        col_idx += 1

    # Build table rows
    rows_html = []
    for row in data:
        cells = []
        isolate = row.get('isolate', '')
        treatment_rec = row.get('treatment recommendation', '')

        cells.append(f'            <td colname="isolate">{isolate}</td>')
        cells.append(f'            <td colname="treatment recommendation">{treatment_rec}</td>')

        for abx in antibiotics:
            nwt_col = f'{abx}_NWT'
            wt_col = f'{abx}_WT'
            nwt_val = row.get(nwt_col, '')
            wt_val = row.get(wt_col, '')

            cells.append(f'            <td colname="{abx}R">{nwt_val}</td>')
            cells.append(f'            <td colname="{abx}S">{wt_val}</td>')

        rows_html.append('          <tr>\n' + '\n'.join(cells) + '\n          </tr>')

    table_html = f'''    <div class="table-wrapper">
      <table id="resultsTable">
        <thead>
          <tr>
{chr(10).join(headers_html)}
          </tr>
        </thead>
        <tbody>
{chr(10).join(rows_html)}
        </tbody>
      </table>
    </div>

    <!-- Pagination Controls -->
    <div class="pagination-controls">
      <div class="page-size-selector">
        <label for="page-size-resultsTable">Show:</label>
        <select id="page-size-resultsTable" onchange="changePageSize(this.value)">
          <option value="10" selected>10</option>
          <option value="25">25</option>
          <option value="50">50</option>
          <option value="100">100</option>
          <option value="500">500</option>
          <option value="1000">1000</option>
          <option value="-1">All</option>
        </select>
      </div>
      <div class="pagination-info" id="pagination-info-resultsTable"></div>
      <div class="pagination-buttons" id="pagination-buttons-resultsTable">
        <button id="prev-page-resultsTable" onclick="previousPage()">Previous</button>
        <button id="next-page-resultsTable" onclick="nextPage()">Next</button>
      </div>
    </div>'''

    return table_html


def _generate_treatment_table_html(data: List[Dict[str, str]]) -> str:
    """
    Generates the treatment recommendations table HTML.

    Args:
        data: List of dictionaries from treatment_output.tsv

    Returns:
        HTML string for treatment table
    """
    if not data:
        return '<p>No treatment data available</p>'

    # Define columns to display (simplified to only 4 essential columns)
    # Excluded: recommended_1, recommended_2, chosen_regimen (internal use)
    # Excluded: individual regimen columns (ceftriaxone+azithromycin, ceftriaxone, etc.)
    columns = [
        'isolate',
        'Predicted Profile',
        'Recommended Treatment',
        'Comment'
    ]

    # Build table headers
    headers_html = []
    for idx, col in enumerate(columns):
        display_name = col.replace('_', ' ').title()
        headers_html.append(f'            <th colname="{col}" onclick="sortTable({idx})">{display_name}</th>')

    # Build table rows
    rows_html = []

    for row in data:
        cells = []
        for col in columns:
            value = row.get(col, '')
            cells.append(f'            <td colname="{col}">{value}</td>')

        rows_html.append('          <tr>\n' + '\n'.join(cells) + '\n          </tr>')

    table_html = f'''    <div class="table-wrapper">
      <table id="treatmentTable">
        <thead>
          <tr>
{chr(10).join(headers_html)}
          </tr>
        </thead>
        <tbody>
{chr(10).join(rows_html)}
        </tbody>
      </table>
    </div>

    <!-- Pagination Controls -->
    <div class="pagination-controls">
      <div class="page-size-selector">
        <label for="page-size-treatmentTable">Show:</label>
        <select id="page-size-treatmentTable" onchange="changePageSize(this.value)">
          <option value="10" selected>10</option>
          <option value="25">25</option>
          <option value="50">50</option>
          <option value="100">100</option>
          <option value="500">500</option>
          <option value="1000">1000</option>
          <option value="-1">All</option>
        </select>
      </div>
      <div class="pagination-info" id="pagination-info-treatmentTable"></div>
      <div class="pagination-buttons" id="pagination-buttons-treatmentTable">
        <button id="prev-page-treatmentTable" onclick="previousPage()">Previous</button>
        <button id="next-page-treatmentTable" onclick="nextPage()">Next</button>
      </div>
    </div>'''

    return table_html


def generate_resistance_profile_html(tsv_path: str, output_html_path: str, antibiotics: List[str] = None) -> None:
    """
    Generates standalone resistance profile HTML (for sensiscript standalone use).

    Args:
        tsv_path: Path to sensiscript TSV output
        output_html_path: Path for output HTML file
        antibiotics: List of antibiotic names (optional, will auto-detect from TSV if not provided)
    """
    # Collapse redundant separators (e.g. "outdir//results.html") so the written and reported paths match
    output_html_path = os.path.normpath(output_html_path)

    # Auto-detect antibiotics from TSV if not provided
    if antibiotics is None or len(antibiotics) == 0:
        antibiotics = _extract_antibiotics_from_tsv(tsv_path)

    # Read TSV data
    data = _read_tsv_to_dict(tsv_path)
    isolate_count = len(data)

    # Generate table HTML
    table_html = _generate_resistance_table_html(data, antibiotics)

    # Generate controls HTML
    antibiotic_options = []
    for abx in antibiotics:
        antibiotic_options.append(f'          <option value="{abx}">{abx.capitalize()}</option>')

    controls_html = f'''    <div class="controls">
      <div class="control-group">
        <label for="show-cols-abx">Filter Columns</label>
        <select id="show-cols-abx" multiple size="4" onchange="filterColumnsByAbx()">
{chr(10).join(antibiotic_options)}
          <option value="all" selected>All antibiotics</option>
        </select>
      </div>

      <div class="control-group" style="flex: 1;">
        <label for="search-box-resultsTable">Search</label>
        <input type="text" id="search-box-resultsTable" onkeyup="searchTable()" placeholder="Filter by mutation, gene, or isolate...">
      </div>
    </div>'''

    # Build complete HTML
    css = _get_css_template()
    js = _get_javascript_template('resultsTable', antibiotics)

    html_content = f'''<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>Sensiscript Results</title>
{css}

{js}
</head>
<body>
  <div class="container">
    <h2><em>Neisseria gonorrhoeae</em> sensityper results</h2>
    <div class="subtitle">Treatment recommendation from genome-based antimicrobial resistance profiles <span class="info-badge">{isolate_count} isolates</span></div>

{controls_html}

    <div class="toggle-switch-container">
      <label class="toggle-switch">
        <span>Show WT columns</span>
        <input type="checkbox" id="show-wt-cols" onchange="toggleWTColumns()" checked>
        <span class="toggle-slider"></span>
      </label>
    </div>

{table_html}
  </div>
</body>
</html>
'''

    # Write to file
    with open(output_html_path, 'w', encoding='utf-8') as f:
        f.write(html_content)

    print(f"Generated resistance profile HTML: {output_html_path}")


def generate_combined_tabbed_html(
    sensiscript_tsv_path: str,
    treatment_tsv_path: str,
    output_html_path: str,
    antibiotics: List[str] = None,
    alert_tsv_path: Optional[str] = None
) -> None:
    """
    Generates combined HTML with tabbed interface (for sensityper).

    Args:
        sensiscript_tsv_path: Path to sensiscript TSV output
        treatment_tsv_path: Path to treatment_output.tsv
        output_html_path: Path for output HTML file
        antibiotics: List of antibiotic names (optional, will auto-detect from TSV if not provided)
        alert_tsv_path: Optional path to alert_output.tsv
    """
    # Collapse redundant separators (e.g. "outdir//treatment.html") so the written and reported paths match
    output_html_path = os.path.normpath(output_html_path)

    # Auto-detect antibiotics from TSV if not provided
    if antibiotics is None or len(antibiotics) == 0:
        antibiotics = _extract_antibiotics_from_tsv(sensiscript_tsv_path)

    # Read data
    resistance_data = _read_tsv_to_dict(sensiscript_tsv_path)
    treatment_data = _read_tsv_to_dict(treatment_tsv_path)
    alert_data = _check_alert_content(alert_tsv_path)

    isolate_count = len(resistance_data)

    # Antibiotics available in this setting (from --available_antibiotics); these gate which
    # regimens sensitreat may choose, so report them alongside the predictions.
    available_html = ''
    if antibiotics:
        chips = ''.join('<span class="info-badge">{a}</span>'.format(a=a) for a in antibiotics)
        available_html = ('\n      <div style="margin-top: 6px;">Antibiotics available for treatment '
                          'in this setting{c}</div>'.format(c=chips))

    # Generate alert banner (always shown, displays "No clinical alerts" if none exist)
    alert_html = _generate_alert_banner_html(alert_data)

    # Generate resistance table
    resistance_table = _generate_resistance_table_html(resistance_data, antibiotics)

    # Generate resistance controls
    antibiotic_options = []
    for abx in antibiotics:
        antibiotic_options.append(f'          <option value="{abx}">{abx.capitalize()}</option>')

    resistance_controls = f'''      <div class="controls">
        <div class="control-group">
          <label for="show-cols-abx">Filter Columns</label>
          <select id="show-cols-abx" multiple size="4" onchange="filterColumnsByAbx()">
{chr(10).join(antibiotic_options)}
            <option value="all" selected>All antibiotics</option>
          </select>
        </div>

        <div class="control-group" style="flex: 1;">
          <label for="search-box-resultsTable">Search</label>
          <input type="text" id="search-box-resultsTable" onkeyup="searchTable()" placeholder="Filter by mutation, gene, or isolate...">
        </div>
      </div>'''

    # Generate treatment table
    treatment_table = _generate_treatment_table_html(treatment_data)

    # Treatment controls (search only)
    treatment_controls = '''      <div class="controls">
        <div class="control-group" style="flex: 1;">
          <label for="search-box-treatmentTable">Search</label>
          <input type="text" id="search-box-treatmentTable" onkeyup="searchTable()" placeholder="Filter by isolate, treatment, or comment...">
        </div>
      </div>'''

    # Build complete HTML
    css = _get_css_template()
    # Use 'treatmentTable' as default table ID for combined HTML (Treatment tab shown first)
    js_combined = _get_javascript_template('treatmentTable', antibiotics)

    html_content = f'''<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>Sensityper Results - Neisseria gonorrhoeae</title>
{css}

{js_combined}
</head>
<body>
  <div class="container">
    <h2><em>Neisseria gonorrhoeae</em> sensityper results</h2>
    <div class="subtitle">Treatment recommendation from genome-based antimicrobial resistance profiles <span class="info-badge">{isolate_count} isolates</span>{available_html}</div>

{alert_html}
    <!-- Tab Navigation -->
    <div class="tab-container">
      <nav class="tab-nav">
        <button class="active" onclick="switchTab('treatment')">Treatment Recommendations</button>
        <button onclick="switchTab('resistance')">Resistance Profile</button>
      </nav>
    </div>

    <!-- Tab 1: Treatment Recommendations (shown first) -->
    <div id="treatment-tab" class="tab-content active">
{treatment_controls}

{treatment_table}
    </div>

    <!-- Tab 2: Resistance Profile -->
    <div id="resistance-tab" class="tab-content">
{resistance_controls}

      <div class="toggle-switch-container">
        <label class="toggle-switch">
          <span>Show WT columns</span>
          <input type="checkbox" id="show-wt-cols" onchange="toggleWTColumns()" checked>
          <span class="toggle-slider"></span>
        </label>
      </div>

{resistance_table}
    </div>
  </div>
</body>
</html>
'''

    # Write to file
    with open(output_html_path, 'w', encoding='utf-8') as f:
        f.write(html_content)

    print(f"Generated combined tabbed HTML: {output_html_path}")

