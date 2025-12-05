"""
Pathogen Report Analyzer - A Flask-based web application for analyzing pathogen detection data
Author: VetLab
Description: Interactive tool for selecting, filtering, and visualizing pathogen detection reports
"""

import os
import sys
from pathlib import Path
from typing import List, Dict, Optional, Tuple
from dataclasses import dataclass
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
from flask import Flask, render_template_string, request, jsonify
from werkzeug.utils import secure_filename
import json
import traceback


# ============================================================================
# DATA MODELS
# ============================================================================

@dataclass
class PathogenData:
    """Data model for pathogen information"""
    pathogen_type: str
    full_pathogen_name: str
    binomial_name: str
    full_taxonomic_tree: str
    similarity_value: float
    used_tools: str
    coverage_ratio: float
    length_of_read: int
    quality_of_read: float
    database: str
    number_of_hits: int
    row_index: int


@dataclass
class ReportFile:
    """Data model for report file information"""
    filename: str
    filepath: str
    sheet_names: List[str]
    has_distilled_report: bool
    file_size: int
    modification_time: str


# ============================================================================
# CORE BUSINESS LOGIC CLASSES
# ============================================================================

class ExcelReportReader:
    """Handles reading and parsing Excel report files"""
    
    REQUIRED_SHEET = "Distilled Report"
    REQUIRED_COLUMNS = [
        'Pathogen Type', 'full pathogen name', 'Binomial name',
        'Full Taxonomic Tree', 'similarity value', 'used tools',
        'coverage ratio', 'Length of that exact read', 'quality of read',
        'database', 'Number of Hits'
    ]
    
    def __init__(self, filepath: str):
        self.filepath = filepath
        self.filename = os.path.basename(filepath)
        
    def validate_file(self) -> Tuple[bool, str]:
        """Validate if file exists and is readable"""
        if not os.path.exists(self.filepath):
            return False, "File does not exist"
        
        if not self.filepath.endswith(('.xlsx', '.xls')):
            return False, "File is not an Excel file"
        
        try:
            pd.ExcelFile(self.filepath)
            return True, "Valid"
        except Exception as e:
            return False, f"Cannot read file: {str(e)}"
    
    def get_sheet_names(self) -> List[str]:
        """Get all sheet names from the Excel file"""
        try:
            excel_file = pd.ExcelFile(self.filepath)
            return excel_file.sheet_names
        except Exception as e:
            print(f"Error reading sheets: {e}")
            return []
    
    def has_distilled_report(self) -> bool:
        """Check if file contains 'Distilled Report' sheet"""
        return self.REQUIRED_SHEET in self.get_sheet_names()
    
    def read_distilled_report(self) -> Optional[pd.DataFrame]:
        """Read the 'Distilled Report' sheet"""
        try:
            df = pd.read_excel(self.filepath, sheet_name=self.REQUIRED_SHEET)
            return df
        except Exception as e:
            print(f"Error reading Distilled Report: {e}")
            return None
    
    def validate_columns(self, df: pd.DataFrame) -> Tuple[bool, List[str]]:
        """Validate if all required columns exist"""
        missing_columns = [col for col in self.REQUIRED_COLUMNS if col not in df.columns]
        return len(missing_columns) == 0, missing_columns
    
    def parse_data(self, df: pd.DataFrame, selected_rows: List[int] = None) -> List[PathogenData]:
        """Parse DataFrame into PathogenData objects"""
        data_objects = []
        
        if selected_rows:
            df = df.iloc[selected_rows]
        
        for idx, row in df.iterrows():
            try:
                data_obj = PathogenData(
                    pathogen_type=str(row.get('Pathogen Type', '')),
                    full_pathogen_name=str(row.get('full pathogen name', '')),
                    binomial_name=str(row.get('Binomial name', '')),
                    full_taxonomic_tree=str(row.get('Full Taxonomic Tree', '')),
                    similarity_value=float(row.get('similarity value', 0)),
                    used_tools=str(row.get('used tools', '')),
                    coverage_ratio=float(row.get('coverage ratio', 0)),
                    length_of_read=int(row.get('Length of that exact read', 0)),
                    quality_of_read=float(row.get('quality of read', 0)),
                    database=str(row.get('database', '')),
                    number_of_hits=int(row.get('Number of Hits', 0)),
                    row_index=int(idx)
                )
                data_objects.append(data_obj)
            except Exception as e:
                print(f"Error parsing row {idx}: {e}")
                continue
        
        return data_objects


class FolderScanner:
    """Scans folders for Excel files"""
    
    EXCEL_EXTENSIONS = ('.xlsx', '.xls')
    
    def __init__(self, folder_path: str):
        self.folder_path = Path(folder_path)
    
    def scan(self) -> List[ReportFile]:
        """Scan folder for Excel files"""
        report_files = []
        
        try:
            for file_path in self.folder_path.glob('*'):
                if file_path.suffix.lower() in self.EXCEL_EXTENSIONS:
                    reader = ExcelReportReader(str(file_path))
                    
                    report_file = ReportFile(
                        filename=file_path.name,
                        filepath=str(file_path),
                        sheet_names=reader.get_sheet_names(),
                        has_distilled_report=reader.has_distilled_report(),
                        file_size=file_path.stat().st_size,
                        modification_time=pd.Timestamp(file_path.stat().st_mtime, unit='s').strftime('%Y-%m-%d %H:%M:%S')
                    )
                    report_files.append(report_file)
        except Exception as e:
            print(f"Error scanning folder: {e}")
        
        return sorted(report_files, key=lambda x: x.filename)


class DataFilterer:
    """Filters pathogen data based on user criteria"""
    
    @staticmethod
    def filter_by_coverage_ratio(data: List[PathogenData], min_ratio: float) -> List[PathogenData]:
        """Filter data by minimum coverage ratio"""
        return [d for d in data if d.coverage_ratio >= min_ratio]
    
    @staticmethod
    def filter_by_number_of_hits(data: List[PathogenData], min_hits: int) -> List[PathogenData]:
        """Filter data by minimum number of hits"""
        return [d for d in data if d.number_of_hits >= min_hits]
    
    @staticmethod
    def filter_by_row_selection(data: List[PathogenData], selected_row_index: int) -> List[PathogenData]:
        """Filter to keep only rows above (and including) the selected row"""
        return [d for d in data if d.row_index <= selected_row_index]


class ChartGenerator:
    """Generates interactive charts using Plotly"""
    
    # Color palettes for attractive visualizations
    COLORS = {
        'vibrant': ['#FF6B6B', '#4ECDC4', '#45B7D1', '#FFA07A', '#98D8C8', '#F7DC6F', '#BB8FCE', '#85C1E2', '#F8B739', '#52B788'],
        'pastel': ['#FFB6C1', '#B0E0E6', '#98FB98', '#DDA0DD', '#F0E68C', '#FFE4B5', '#E0BBE4', '#FFDAB9', '#D4F1F4', '#C5E1A5'],
        'professional': ['#2E86AB', '#A23B72', '#F18F01', '#C73E1D', '#6A994E', '#BC4B51', '#5E548E', '#277DA1', '#F77F00', '#06A77D']
    }
    
    def __init__(self, color_scheme: str = 'vibrant'):
        self.colors = self.COLORS.get(color_scheme, self.COLORS['vibrant'])
    
    def create_pie_chart_coverage(self, data: List[PathogenData]) -> str:
        """Create pie chart for coverage ratio distribution by pathogen"""
        if not data:
            return self._create_empty_chart("No data available for pie chart")
        
        # Aggregate coverage by pathogen type
        pathogen_coverage = {}
        for item in data:
            pathogen = item.full_pathogen_name or item.pathogen_type
            if pathogen in pathogen_coverage:
                pathogen_coverage[pathogen] += item.coverage_ratio
            else:
                pathogen_coverage[pathogen] = item.coverage_ratio
        
        labels = list(pathogen_coverage.keys())
        values = list(pathogen_coverage.values())
        
        fig = go.Figure(data=[go.Pie(
            labels=labels,
            values=values,
            hole=0.3,
            marker=dict(colors=self.colors[:len(labels)]),
            textinfo='label+percent',
            textposition='outside',
            hovertemplate='<b>%{label}</b><br>Coverage: %{value:.2f}<br>Percentage: %{percent}<extra></extra>'
        )])
        
        fig.update_layout(
            title={
                'text': 'Coverage Ratio Distribution by Pathogen',
                'x': 0.5,
                'xanchor': 'center',
                'font': {'size': 20, 'color': '#2c3e50', 'family': 'Arial Black'}
            },
            showlegend=True,
            height=600,
            legend=dict(
                orientation="v",
                yanchor="middle",
                y=0.5,
                xanchor="left",
                x=1.1,
                bgcolor='rgba(255,255,255,0.8)',
                bordercolor='#E0E0E0',
                borderwidth=2
            ),
            paper_bgcolor='#F8F9FA',
            plot_bgcolor='#F8F9FA'
        )
        
        return fig.to_html(full_html=False, include_plotlyjs='cdn')
    
    def create_pie_chart_hits(self, data: List[PathogenData]) -> str:
        """Create pie chart for number of hits distribution"""
        if not data:
            return self._create_empty_chart("No data available for pie chart")
        
        # Aggregate hits by pathogen type
        pathogen_hits = {}
        for item in data:
            pathogen = item.full_pathogen_name or item.pathogen_type
            if pathogen in pathogen_hits:
                pathogen_hits[pathogen] += item.number_of_hits
            else:
                pathogen_hits[pathogen] = item.number_of_hits
        
        labels = list(pathogen_hits.keys())
        values = list(pathogen_hits.values())
        
        fig = go.Figure(data=[go.Pie(
            labels=labels,
            values=values,
            hole=0.3,
            marker=dict(colors=self.colors[:len(labels)]),
            textinfo='label+percent',
            textposition='outside',
            hovertemplate='<b>%{label}</b><br>Hits: %{value}<br>Percentage: %{percent}<extra></extra>'
        )])
        
        fig.update_layout(
            title={
                'text': 'Number of Hits Distribution by Pathogen',
                'x': 0.5,
                'xanchor': 'center',
                'font': {'size': 20, 'color': '#2c3e50', 'family': 'Arial Black'}
            },
            showlegend=True,
            height=600,
            legend=dict(
                orientation="v",
                yanchor="middle",
                y=0.5,
                xanchor="left",
                x=1.1,
                bgcolor='rgba(255,255,255,0.8)',
                bordercolor='#E0E0E0',
                borderwidth=2
            ),
            paper_bgcolor='#F8F9FA',
            plot_bgcolor='#F8F9FA'
        )
        
        return fig.to_html(full_html=False, include_plotlyjs='cdn')
    
    def create_stacked_bar_chart(self, data: List[PathogenData]) -> str:
        """Create stacked bar chart showing multiple metrics"""
        if not data:
            return self._create_empty_chart("No data available for bar chart")
        
        # Prepare data for stacked bar chart
        pathogens = [item.full_pathogen_name or item.pathogen_type for item in data]
        coverage_ratios = [item.coverage_ratio for item in data]
        hit_counts = [item.number_of_hits for item in data]
        quality_scores = [item.quality_of_read for item in data]
        similarity_values = [item.similarity_value for item in data]
        
        fig = go.Figure()
        
        # Add traces for each metric
        fig.add_trace(go.Bar(
            name='Coverage Ratio',
            x=pathogens,
            y=coverage_ratios,
            marker_color=self.colors[0],
            hovertemplate='<b>%{x}</b><br>Coverage Ratio: %{y:.2f}<extra></extra>'
        ))
        
        fig.add_trace(go.Bar(
            name='Number of Hits (scaled)',
            x=pathogens,
            y=[h/10 for h in hit_counts],  # Scale down for visualization
            marker_color=self.colors[1],
            hovertemplate='<b>%{x}</b><br>Number of Hits: %{y:.0f}<extra></extra>'
        ))
        
        fig.add_trace(go.Bar(
            name='Quality Score',
            x=pathogens,
            y=quality_scores,
            marker_color=self.colors[2],
            hovertemplate='<b>%{x}</b><br>Quality Score: %{y:.2f}<extra></extra>'
        ))
        
        fig.add_trace(go.Bar(
            name='Similarity Value',
            x=pathogens,
            y=similarity_values,
            marker_color=self.colors[3],
            hovertemplate='<b>%{x}</b><br>Similarity: %{y:.2f}<extra></extra>'
        ))
        
        fig.update_layout(
            title={
                'text': 'Comprehensive Pathogen Metrics Analysis',
                'x': 0.5,
                'xanchor': 'center',
                'font': {'size': 20, 'color': '#2c3e50', 'family': 'Arial Black'}
            },
            barmode='stack',
            xaxis_title='Pathogen',
            yaxis_title='Metric Value',
            height=600,
            showlegend=True,
            legend=dict(
                orientation="h",
                yanchor="bottom",
                y=1.02,
                xanchor="right",
                x=1,
                bgcolor='rgba(255,255,255,0.8)',
                bordercolor='#E0E0E0',
                borderwidth=2
            ),
            paper_bgcolor='#F8F9FA',
            plot_bgcolor='#FFFFFF',
            xaxis=dict(
                tickangle=-45,
                gridcolor='#E0E0E0'
            ),
            yaxis=dict(
                gridcolor='#E0E0E0'
            )
        )
        
        return fig.to_html(full_html=False, include_plotlyjs='cdn')
    
    def create_quality_heatmap(self, data: List[PathogenData]) -> str:
        """Create heatmap showing quality metrics across pathogens"""
        if not data:
            return self._create_empty_chart("No data available for heatmap")
        
        # Prepare data matrix
        pathogens = [item.full_pathogen_name or item.pathogen_type for item in data]
        metrics = ['Coverage Ratio', 'Quality Score', 'Similarity Value', 'Read Length (scaled)']
        
        data_matrix = [
            [item.coverage_ratio for item in data],
            [item.quality_of_read for item in data],
            [item.similarity_value for item in data],
            [item.length_of_read / 1000 for item in data]  # Scale down
        ]
        
        fig = go.Figure(data=go.Heatmap(
            z=data_matrix,
            x=pathogens,
            y=metrics,
            colorscale='Viridis',
            hovertemplate='<b>%{y}</b><br>Pathogen: %{x}<br>Value: %{z:.2f}<extra></extra>'
        ))
        
        fig.update_layout(
            title={
                'text': 'Quality Metrics Heatmap',
                'x': 0.5,
                'xanchor': 'center',
                'font': {'size': 20, 'color': '#2c3e50', 'family': 'Arial Black'}
            },
            xaxis_title='Pathogen',
            yaxis_title='Metric',
            height=500,
            paper_bgcolor='#F8F9FA',
            plot_bgcolor='#F8F9FA'
        )
        
        return fig.to_html(full_html=False, include_plotlyjs='cdn')
    
    def create_scatter_plot(self, data: List[PathogenData]) -> str:
        """Create scatter plot of coverage vs hits"""
        if not data:
            return self._create_empty_chart("No data available for scatter plot")
        
        pathogens = [item.full_pathogen_name or item.pathogen_type for item in data]
        coverage = [item.coverage_ratio for item in data]
        hits = [item.number_of_hits for item in data]
        quality = [item.quality_of_read for item in data]
        
        fig = go.Figure()
        
        fig.add_trace(go.Scatter(
            x=coverage,
            y=hits,
            mode='markers',
            marker=dict(
                size=[q*2 for q in quality],  # Size based on quality
                color=quality,
                colorscale='Plasma',
                showscale=True,
                colorbar=dict(title="Quality Score"),
                line=dict(width=1, color='white')
            ),
            text=pathogens,
            hovertemplate='<b>%{text}</b><br>Coverage: %{x:.2f}<br>Hits: %{y}<br>Quality: %{marker.color:.2f}<extra></extra>'
        ))
        
        fig.update_layout(
            title={
                'text': 'Coverage Ratio vs Number of Hits',
                'x': 0.5,
                'xanchor': 'center',
                'font': {'size': 20, 'color': '#2c3e50', 'family': 'Arial Black'}
            },
            xaxis_title='Coverage Ratio',
            yaxis_title='Number of Hits',
            height=600,
            paper_bgcolor='#F8F9FA',
            plot_bgcolor='#FFFFFF',
            xaxis=dict(gridcolor='#E0E0E0'),
            yaxis=dict(gridcolor='#E0E0E0')
        )
        
        return fig.to_html(full_html=False, include_plotlyjs='cdn')
    
    def _create_empty_chart(self, message: str) -> str:
        """Create an empty chart with a message"""
        fig = go.Figure()
        fig.add_annotation(
            text=message,
            xref="paper",
            yref="paper",
            x=0.5,
            y=0.5,
            showarrow=False,
            font=dict(size=16, color="gray")
        )
        fig.update_layout(
            xaxis=dict(showgrid=False, showticklabels=False, zeroline=False),
            yaxis=dict(showgrid=False, showticklabels=False, zeroline=False),
            height=400
        )
        return fig.to_html(full_html=False, include_plotlyjs='cdn')


# ============================================================================
# FLASK APPLICATION
# ============================================================================

class PathogenReportAnalyzerApp:
    """Main application class"""
    
    def __init__(self):
        self.app = Flask(__name__)
        self.app.config['MAX_CONTENT_LENGTH'] = 100 * 1024 * 1024  # 100MB max
        self.current_folder = None
        self.current_file = None
        self.current_data = None
        self.chart_generator = ChartGenerator(color_scheme='vibrant')
        self._setup_routes()
    
    def _setup_routes(self):
        """Setup Flask routes"""
        self.app.route('/')(self.index)
        self.app.route('/browse_folders', methods=['GET'])(self.browse_folders)
        self.app.route('/select_folder', methods=['POST'])(self.select_folder)
        self.app.route('/scan_files', methods=['GET'])(self.scan_files)
        self.app.route('/load_file', methods=['POST'])(self.load_file)
        self.app.route('/get_file_data', methods=['GET'])(self.get_file_data)
        self.app.route('/filter_and_plot', methods=['POST'])(self.filter_and_plot)
        self.app.route('/health', methods=['GET'])(self.health)
    
    def index(self):
        """Main page"""
        return render_template_string(HTML_TEMPLATE)
    
    def browse_folders(self):
        """Browse folders endpoint"""
        path = request.args.get('path', os.path.expanduser('~'))
        
        try:
            path = os.path.expanduser(path)
            path = os.path.abspath(path)
            
            folders = []
            parent = str(Path(path).parent)
            
            if path != parent:
                folders.append({
                    'name': '.. (Parent Directory)',
                    'path': parent,
                    'is_parent': True
                })
            
            try:
                items = sorted(os.listdir(path))
            except PermissionError:
                return jsonify({'error': 'Permission denied'}), 403
            
            for item in items:
                item_path = os.path.join(path, item)
                try:
                    if not item.startswith('.') and os.path.isdir(item_path):
                        folders.append({
                            'name': item,
                            'path': item_path,
                            'is_parent': False
                        })
                except (PermissionError, OSError):
                    continue
            
            return jsonify({'folders': folders, 'current_path': path})
        except Exception as e:
            return jsonify({'error': str(e)}), 400
    
    def select_folder(self):
        """Select folder endpoint"""
        try:
            folder_path = request.json.get('path')
            
            if not folder_path:
                return jsonify({'success': False, 'error': 'No path provided'}), 400
            
            folder_path = os.path.expanduser(folder_path)
            folder_path = os.path.abspath(folder_path)
            
            if os.path.isdir(folder_path):
                self.current_folder = folder_path
                print(f"\n✓ Selected folder: {self.current_folder}\n")
                return jsonify({'success': True, 'path': self.current_folder})
            else:
                return jsonify({'success': False, 'error': 'Invalid directory'}), 400
        except Exception as e:
            return jsonify({'success': False, 'error': str(e)}), 400
    
    def scan_files(self):
        """Scan files in selected folder"""
        if not self.current_folder:
            return jsonify({'error': 'No folder selected'}), 400
        
        try:
            scanner = FolderScanner(self.current_folder)
            files = scanner.scan()
            
            files_data = [
                {
                    'filename': f.filename,
                    'filepath': f.filepath,
                    'has_distilled_report': f.has_distilled_report,
                    'file_size': f'{f.file_size / 1024:.2f} KB',
                    'modification_time': f.modification_time,
                    'sheet_count': len(f.sheet_names)
                }
                for f in files
            ]
            
            return jsonify({
                'success': True,
                'files': files_data,
                'count': len(files_data)
            })
        except Exception as e:
            return jsonify({'error': str(e)}), 400
    
    def load_file(self):
        """Load selected Excel file"""
        try:
            filepath = request.json.get('filepath')
            
            if not filepath:
                return jsonify({'success': False, 'error': 'No file specified'}), 400
            
            reader = ExcelReportReader(filepath)
            
            # Validate file
            is_valid, message = reader.validate_file()
            if not is_valid:
                return jsonify({'success': False, 'error': message}), 400
            
            # Read data
            df = reader.read_distilled_report()
            if df is None:
                return jsonify({'success': False, 'error': 'Could not read Distilled Report sheet'}), 400
            
            # Validate columns
            is_valid, missing_cols = reader.validate_columns(df)
            if not is_valid:
                return jsonify({
                    'success': False,
                    'error': f'Missing columns: {", ".join(missing_cols)}'
                }), 400
            
            # Parse data
            self.current_data = reader.parse_data(df)
            self.current_file = filepath
            
            # Prepare summary
            summary = {
                'filename': os.path.basename(filepath),
                'total_rows': len(self.current_data),
                'columns': list(df.columns),
                'preview_data': df.head(10).to_dict('records'),
                'coverage_range': {
                    'min': df['coverage ratio'].min(),
                    'max': df['coverage ratio'].max(),
                    'mean': df['coverage ratio'].mean()
                },
                'hits_range': {
                    'min': int(df['Number of Hits'].min()),
                    'max': int(df['Number of Hits'].max()),
                    'mean': float(df['Number of Hits'].mean())
                }
            }
            
            return jsonify({
                'success': True,
                'summary': summary
            })
            
        except Exception as e:
            traceback.print_exc()
            return jsonify({'success': False, 'error': str(e)}), 400
    
    def get_file_data(self):
        """Get current file data"""
        if not self.current_data:
            return jsonify({'error': 'No file loaded'}), 400
        
        try:
            data_dicts = [
                {
                    'row_index': d.row_index,
                    'pathogen_type': d.pathogen_type,
                    'full_pathogen_name': d.full_pathogen_name,
                    'binomial_name': d.binomial_name,
                    'coverage_ratio': d.coverage_ratio,
                    'number_of_hits': d.number_of_hits,
                    'quality_of_read': d.quality_of_read,
                    'similarity_value': d.similarity_value
                }
                for d in self.current_data
            ]
            
            return jsonify({
                'success': True,
                'data': data_dicts
            })
        except Exception as e:
            return jsonify({'error': str(e)}), 400
    
    def filter_and_plot(self):
        """Filter data and generate plots"""
        if not self.current_data:
            return jsonify({'error': 'No file loaded'}), 400
        
        try:
            selected_row_index = request.json.get('selected_row_index')
            
            if selected_row_index is None:
                return jsonify({'error': 'No row selected'}), 400
            
            # Filter data (keep rows from 0 to selected_row_index inclusive)
            filtered_data = DataFilterer.filter_by_row_selection(
                self.current_data,
                selected_row_index
            )
            
            if not filtered_data:
                return jsonify({'error': 'No data after filtering'}), 400
            
            # Generate charts
            pie_chart_coverage = self.chart_generator.create_pie_chart_coverage(filtered_data)
            pie_chart_hits = self.chart_generator.create_pie_chart_hits(filtered_data)
            stacked_bar = self.chart_generator.create_stacked_bar_chart(filtered_data)
            heatmap = self.chart_generator.create_quality_heatmap(filtered_data)
            scatter = self.chart_generator.create_scatter_plot(filtered_data)
            
            return jsonify({
                'success': True,
                'charts': {
                    'pie_coverage': pie_chart_coverage,
                    'pie_hits': pie_chart_hits,
                    'stacked_bar': stacked_bar,
                    'heatmap': heatmap,
                    'scatter': scatter
                },
                'filtered_count': len(filtered_data)
            })
            
        except Exception as e:
            traceback.print_exc()
            return jsonify({'error': str(e)}), 400
    
    def health(self):
        """Health check"""
        return jsonify({
            'status': 'running',
            'current_folder': self.current_folder,
            'current_file': self.current_file,
            'data_loaded': self.current_data is not None
        })
    
    def run(self, host='0.0.0.0', port=5000, debug=True):
        """Run the Flask application"""
        self.app.run(host=host, port=port, debug=debug, use_reloader=False)


# ============================================================================
# HTML TEMPLATE
# ============================================================================

HTML_TEMPLATE = """
<!DOCTYPE html>
<html>
<head>
    <title>Pathogen Report Analyzer</title>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
    <style>
        * {
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }
        
        body {
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            min-height: 100vh;
            padding: 20px;
        }
        
        .container {
            max-width: 1400px;
            margin: 0 auto;
            background: white;
            border-radius: 15px;
            padding: 30px;
            box-shadow: 0 20px 60px rgba(0,0,0,0.3);
        }
        
        h1 {
            text-align: center;
            color: #2c3e50;
            font-size: 2.5em;
            margin-bottom: 10px;
            text-shadow: 2px 2px 4px rgba(0,0,0,0.1);
        }
        
        .subtitle {
            text-align: center;
            color: #7f8c8d;
            margin-bottom: 30px;
            font-size: 1.1em;
        }
        
        .step {
            background: #f8f9fa;
            border-left: 5px solid #4CAF50;
            padding: 20px;
            margin-bottom: 20px;
            border-radius: 8px;
        }
        
        .step-title {
            font-size: 1.3em;
            color: #2c3e50;
            margin-bottom: 15px;
            display: flex;
            align-items: center;
            gap: 10px;
        }
        
        .step-number {
            background: #4CAF50;
            color: white;
            width: 35px;
            height: 35px;
            border-radius: 50%;
            display: flex;
            align-items: center;
            justify-content: center;
            font-weight: bold;
        }
        
        .folder-browser {
            display: grid;
            grid-template-columns: 2fr 1fr;
            gap: 20px;
            margin-top: 15px;
        }
        
        .current-path {
            background: white;
            padding: 12px;
            border-radius: 6px;
            border: 2px solid #e0e0e0;
            font-family: monospace;
            font-size: 0.9em;
            word-break: break-all;
        }
        
        .folder-list {
            max-height: 300px;
            overflow-y: auto;
            border: 2px solid #e0e0e0;
            border-radius: 6px;
            background: white;
            margin: 15px 0;
        }
        
        .folder-item, .file-item {
            padding: 12px 15px;
            cursor: pointer;
            border-bottom: 1px solid #e0e0e0;
            transition: all 0.2s;
            display: flex;
            align-items: center;
            gap: 10px;
        }
        
        .folder-item:hover, .file-item:hover {
            background: #f0f7ff;
            padding-left: 20px;
        }
        
        .folder-item.parent {
            background: #e3f2fd;
            font-weight: bold;
        }
        
        .file-item.has-report {
            border-left: 4px solid #4CAF50;
        }
        
        .file-item.selected {
            background: #e8f5e9;
            border-left: 4px solid #2e7d32;
        }
        
        button {
            background: #4CAF50;
            color: white;
            padding: 12px 24px;
            border: none;
            border-radius: 6px;
            cursor: pointer;
            font-size: 1em;
            font-weight: bold;
            transition: all 0.3s;
            box-shadow: 0 2px 5px rgba(0,0,0,0.2);
        }
        
        button:hover {
            background: #45a049;
            transform: translateY(-2px);
            box-shadow: 0 4px 8px rgba(0,0,0,0.3);
        }
        
        button:disabled {
            background: #cccccc;
            cursor: not-allowed;
            transform: none;
        }
        
        .btn-secondary {
            background: #2196F3;
        }
        
        .btn-secondary:hover {
            background: #1976D2;
        }
        
        .data-table {
            width: 100%;
            border-collapse: collapse;
            margin-top: 15px;
            background: white;
            border-radius: 8px;
            overflow: hidden;
            box-shadow: 0 2px 8px rgba(0,0,0,0.1);
        }
        
        .data-table th {
            background: #4CAF50;
            color: white;
            padding: 12px;
            text-align: left;
            font-weight: bold;
        }
        
        .data-table td {
            padding: 10px;
            border-bottom: 1px solid #e0e0e0;
        }
        
        .data-table tr:hover {
            background: #f5f5f5;
            cursor: pointer;
        }
        
        .data-table tr.selected {
            background: #e8f5e9;
            border-left: 4px solid #4CAF50;
        }
        
        .charts-container {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(500px, 1fr));
            gap: 20px;
            margin-top: 20px;
        }
        
        .chart-box {
            background: white;
            border-radius: 8px;
            padding: 15px;
            box-shadow: 0 2px 8px rgba(0,0,0,0.1);
        }
        
        .loading {
            text-align: center;
            padding: 20px;
            color: #666;
        }
        
        .error {
            background: #ffebee;
            color: #c62828;
            padding: 15px;
            border-radius: 6px;
            margin-top: 15px;
            display: none;
            border-left: 4px solid #c62828;
        }
        
        .success {
            background: #e8f5e9;
            color: #2e7d32;
            padding: 15px;
            border-radius: 6px;
            margin-top: 15px;
            display: none;
            border-left: 4px solid #4CAF50;
        }
        
        .info-box {
            background: #e3f2fd;
            padding: 15px;
            border-radius: 6px;
            margin-top: 15px;
            border-left: 4px solid #2196F3;
        }
        
        .badge {
            display: inline-block;
            padding: 4px 8px;
            border-radius: 4px;
            font-size: 0.85em;
            font-weight: bold;
        }
        
        .badge-success {
            background: #4CAF50;
            color: white;
        }
        
        .badge-warning {
            background: #ff9800;
            color: white;
        }
        
        .hidden {
            display: none;
        }
        
        @keyframes fadeIn {
            from { opacity: 0; transform: translateY(20px); }
            to { opacity: 1; transform: translateY(0); }
        }
        
        .fade-in {
            animation: fadeIn 0.5s ease-out;
        }
    </style>
</head>
<body>
    <div class="container">
        <h1>🧬 Pathogen Report Analyzer</h1>
        <p class="subtitle">Advanced Analysis Tool for Pathogen Detection Reports</p>
        
        <!-- STEP 1: Folder Selection -->
        <div class="step" id="step1">
            <div class="step-title">
                <span class="step-number">1</span>
                <span>Select Folder Containing Excel Reports</span>
            </div>
            
            <div class="folder-browser">
                <div>
                    <div class="current-path" id="currentPath">Loading...</div>
                    <div class="folder-list" id="folderList">
                        <div class="loading">Loading folders...</div>
                    </div>
                </div>
                <div>
                    <button onclick="goHome()" style="width: 100%; margin-bottom: 10px;">
                        🏠 Go Home
                    </button>
                    <button onclick="selectFolder()" style="width: 100%;">
                        ✓ Select This Folder
                    </button>
                </div>
            </div>
            
            <div class="error" id="error1"></div>
            <div class="success" id="success1"></div>
        </div>
        
        <!-- STEP 2: File Selection -->
        <div class="step hidden" id="step2">
            <div class="step-title">
                <span class="step-number">2</span>
                <span>Select Excel File to Analyze</span>
            </div>
            
            <button onclick="scanFiles()" class="btn-secondary" style="margin-bottom: 15px;">
                🔍 Scan for Excel Files
            </button>
            
            <div class="folder-list" id="fileList">
                <div class="loading">Click "Scan for Excel Files" to begin</div>
            </div>
            
            <div class="error" id="error2"></div>
            <div class="success" id="success2"></div>
        </div>
        
        <!-- STEP 3: Data Preview and Row Selection -->
        <div class="step hidden" id="step3">
            <div class="step-title">
                <span class="step-number">3</span>
                <span>Select Row Threshold (All rows above will be included)</span>
            </div>
            
            <div class="info-box" id="fileInfo"></div>
            
            <div style="overflow-x: auto; margin-top: 15px;">
                <table class="data-table" id="dataTable">
                    <thead>
                        <tr>
                            <th>Row</th>
                            <th>Pathogen Type</th>
                            <th>Full Name</th>
                            <th>Coverage Ratio</th>
                            <th>Number of Hits</th>
                            <th>Quality</th>
                            <th>Similarity</th>
                        </tr>
                    </thead>
                    <tbody id="dataTableBody">
                    </tbody>
                </table>
            </div>
            
            <button onclick="generatePlots()" style="margin-top: 20px; width: 100%;">
                📊 Generate Plots for Selected Data
            </button>
            
            <div class="error" id="error3"></div>
            <div class="success" id="success3"></div>
        </div>
        
        <!-- STEP 4: Visualization -->
        <div class="step hidden" id="step4">
            <div class="step-title">
                <span class="step-number">4</span>
                <span>Analysis Results & Visualizations</span>
            </div>
            
            <div class="info-box" id="plotInfo"></div>
            
            <div class="charts-container" id="chartsContainer">
            </div>
        </div>
    </div>

    <script>
        let currentPath = '/';
        let homePath = '';
        let selectedFolder = null;
        let selectedFile = null;
        let selectedRowIndex = null;
        let currentData = null;
        
        // Initialize
        window.onload = function() {
            fetch('/health')
                .then(response => response.json())
                .then(data => {
                    homePath = data.current_folder || '/';
                    loadFolders(homePath);
                });
        };
        
        function showError(stepId, message) {
            const errorDiv = document.getElementById(`error${stepId}`);
            errorDiv.textContent = '⚠️ ' + message;
            errorDiv.style.display = 'block';
            setTimeout(() => errorDiv.style.display = 'none', 5000);
        }
        
        function showSuccess(stepId, message) {
            const successDiv = document.getElementById(`success${stepId}`);
            successDiv.textContent = '✓ ' + message;
            successDiv.style.display = 'block';
        }
        
        function goHome() {
            loadFolders(homePath);
        }
        
        function loadFolders(path) {
            const folderList = document.getElementById('folderList');
            folderList.innerHTML = '<div class="loading">Loading folders...</div>';
            
            fetch('/browse_folders?path=' + encodeURIComponent(path))
                .then(response => response.json())
                .then(data => {
                    if (data.error) {
                        showError(1, data.error);
                        return;
                    }
                    
                    currentPath = data.current_path;
                    document.getElementById('currentPath').textContent = currentPath;
                    
                    folderList.innerHTML = '';
                    
                    if (data.folders.length === 0) {
                        folderList.innerHTML = '<div class="loading">No subfolders found</div>';
                        return;
                    }
                    
                    data.folders.forEach(folder => {
                        const div = document.createElement('div');
                        div.className = 'folder-item' + (folder.is_parent ? ' parent' : '');
                        div.innerHTML = `
                            <span>${folder.is_parent ? '⬆️' : '📁'}</span>
                            <span>${folder.name}</span>
                        `;
                        div.onclick = () => loadFolders(folder.path);
                        folderList.appendChild(div);
                    });
                })
                .catch(error => {
                    showError(1, 'Network error: ' + error.message);
                });
        }
        
        function selectFolder() {
            fetch('/select_folder', {
                method: 'POST',
                headers: {'Content-Type': 'application/json'},
                body: JSON.stringify({path: currentPath})
            })
            .then(response => response.json())
            .then(data => {
                if (data.success) {
                    selectedFolder = data.path;
                    showSuccess(1, 'Folder selected: ' + data.path);
                    document.getElementById('step2').classList.remove('hidden');
                    document.getElementById('step2').scrollIntoView({behavior: 'smooth'});
                } else {
                    showError(1, data.error);
                }
            })
            .catch(error => showError(1, 'Error: ' + error.message));
        }
        
        function scanFiles() {
            if (!selectedFolder) {
                showError(2, 'Please select a folder first');
                return;
            }
            
            const fileList = document.getElementById('fileList');
            fileList.innerHTML = '<div class="loading">Scanning for Excel files...</div>';
            
            fetch('/scan_files')
                .then(response => response.json())
                .then(data => {
                    if (data.error) {
                        showError(2, data.error);
                        return;
                    }
                    
                    fileList.innerHTML = '';
                    
                    if (data.files.length === 0) {
                        fileList.innerHTML = '<div class="loading">No Excel files found</div>';
                        return;
                    }
                    
                    data.files.forEach(file => {
                        const div = document.createElement('div');
                        div.className = 'file-item' + (file.has_distilled_report ? ' has-report' : '');
                        div.innerHTML = `
                            <span>📄</span>
                            <div style="flex: 1;">
                                <div><strong>${file.filename}</strong></div>
                                <div style="font-size: 0.85em; color: #666;">
                                    ${file.file_size} | Modified: ${file.modification_time} | 
                                    ${file.has_distilled_report ? 
                                        '<span class="badge badge-success">Has Distilled Report</span>' : 
                                        '<span class="badge badge-warning">No Distilled Report</span>'}
                                </div>
                            </div>
                        `;
                        div.onclick = () => loadFile(file.filepath, div);
                        fileList.appendChild(div);
                    });
                    
                    showSuccess(2, `Found ${data.files.length} Excel file(s)`);
                })
                .catch(error => showError(2, 'Error: ' + error.message));
        }
        
        function loadFile(filepath, element) {
            // Highlight selected file
            document.querySelectorAll('.file-item').forEach(item => {
                item.classList.remove('selected');
            });
            element.classList.add('selected');
            
            selectedFile = filepath;
            
            const step3 = document.getElementById('step3');
            step3.classList.add('hidden');
            
            fetch('/load_file', {
                method: 'POST',
                headers: {'Content-Type': 'application/json'},
                body: JSON.stringify({filepath: filepath})
            })
            .then(response => response.json())
            .then(data => {
                if (data.error) {
                    showError(2, data.error);
                    return;
                }
                
                if (data.success) {
                    showSuccess(2, 'File loaded successfully!');
                    displayFileData(data.summary);
                    step3.classList.remove('hidden');
                    step3.classList.add('fade-in');
                    step3.scrollIntoView({behavior: 'smooth'});
                }
            })
            .catch(error => showError(2, 'Error loading file: ' + error.message));
        }
        
        function displayFileData(summary) {
            const fileInfo = document.getElementById('fileInfo');
            fileInfo.innerHTML = `
                <h3 style="margin-bottom: 10px;">📊 ${summary.filename}</h3>
                <p><strong>Total Rows:</strong> ${summary.total_rows}</p>
                <p><strong>Coverage Ratio Range:</strong> ${summary.coverage_range.min.toFixed(2)} - ${summary.coverage_range.max.toFixed(2)} (Avg: ${summary.coverage_range.mean.toFixed(2)})</p>
                <p><strong>Number of Hits Range:</strong> ${summary.hits_range.min} - ${summary.hits_range.max} (Avg: ${summary.hits_range.mean.toFixed(1)})</p>
            `;
            
            const tableBody = document.getElementById('dataTableBody');
            tableBody.innerHTML = '';
            
            summary.preview_data.forEach((row, index) => {
                const tr = document.createElement('tr');
                tr.innerHTML = `
                    <td>${index}</td>
                    <td>${row['Pathogen Type'] || ''}</td>
                    <td>${row['full pathogen name'] || ''}</td>
                    <td>${typeof row['coverage ratio'] === 'number' ? row['coverage ratio'].toFixed(2) : ''}</td>
                    <td>${row['Number of Hits'] || ''}</td>
                    <td>${typeof row['quality of read'] === 'number' ? row['quality of read'].toFixed(2) : ''}</td>
                    <td>${typeof row['similarity value'] === 'number' ? row['similarity value'].toFixed(2) : ''}</td>
                `;
                tr.onclick = () => selectRow(index, tr);
                tableBody.appendChild(tr);
            });
        }
        
        function selectRow(rowIndex, element) {
            // Highlight selected row
            document.querySelectorAll('.data-table tr').forEach(row => {
                row.classList.remove('selected');
            });
            element.classList.add('selected');
            
            selectedRowIndex = rowIndex;
            showSuccess(3, `Selected row ${rowIndex}. All rows from 0 to ${rowIndex} will be included in analysis.`);
        }
        
        function generatePlots() {
            if (selectedRowIndex === null) {
                showError(3, 'Please select a row threshold first');
                return;
            }
            
            const step4 = document.getElementById('step4');
            step4.classList.add('hidden');
            
            const chartsContainer = document.getElementById('chartsContainer');
            chartsContainer.innerHTML = '<div class="loading">Generating plots...</div>';
            
            fetch('/filter_and_plot', {
                method: 'POST',
                headers: {'Content-Type': 'application/json'},
                body: JSON.stringify({selected_row_index: selectedRowIndex})
            })
            .then(response => response.json())
            .then(data => {
                if (data.error) {
                    showError(3, data.error);
                    return;
                }
                
                if (data.success) {
                    showSuccess(3, `Analysis complete! Generated plots for ${data.filtered_count} rows.`);
                    
                    const plotInfo = document.getElementById('plotInfo');
                    plotInfo.innerHTML = `
                        <h3 style="margin-bottom: 10px;">📈 Analysis Results</h3>
                        <p><strong>Rows Analyzed:</strong> ${data.filtered_count} (from row 0 to ${selectedRowIndex})</p>
                        <p><strong>Charts Generated:</strong> 5 interactive visualizations</p>
                    `;
                    
                    displayCharts(data.charts);
                    step4.classList.remove('hidden');
                    step4.classList.add('fade-in');
                    step4.scrollIntoView({behavior: 'smooth'});
                }
            })
            .catch(error => showError(3, 'Error generating plots: ' + error.message));
        }
        
        function displayCharts(charts) {
            const container = document.getElementById('chartsContainer');
            container.innerHTML = '';
            
            const chartTitles = {
                'pie_coverage': '🥧 Coverage Ratio Distribution',
                'pie_hits': '🎯 Number of Hits Distribution',
                'stacked_bar': '📊 Comprehensive Metrics Analysis',
                'heatmap': '🔥 Quality Metrics Heatmap',
                'scatter': '🔍 Coverage vs Hits Scatter Plot'
            };
            
            Object.entries(charts).forEach(([key, html]) => {
                const chartBox = document.createElement('div');
                chartBox.className = 'chart-box fade-in';
                chartBox.innerHTML = html;
                container.appendChild(chartBox);
            });
        }
    </script>
</body>
</html>
"""


# ============================================================================
# MAIN EXECUTION
# ============================================================================

def main():
    """Main function to run the application"""
    import argparse
    
    parser = argparse.ArgumentParser(
        description='Pathogen Report Analyzer - Advanced Analysis Tool',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python report_plotter.py
  python report_plotter.py --port 8080
  python report_plotter.py --host 127.0.0.1 --port 5000
        """
    )
    
    parser.add_argument(
        '--wdir',
        action='store_true',
        help='Enable working directory selection mode'
    )
    
    parser.add_argument(
        '--port',
        type=int,
        default=5000,
        help='Port to run the server on (default: 5000)'
    )
    
    parser.add_argument(
        '--host',
        type=str,
        default='0.0.0.0',
        help='Host to run the server on (default: 0.0.0.0)'
    )
    
    parser.add_argument(
        '--debug',
        action='store_true',
        help='Run in debug mode'
    )
    
    args = parser.parse_args()
    
    # Print banner
    print("\n" + "="*70)
    print("🧬 PATHOGEN REPORT ANALYZER - Advanced Analysis Tool")
    print("="*70)
    print(f"Version: 2.0.0")
    print(f"Server: http://localhost:{args.port}")
    print(f"Network: http://{args.host}:{args.port}")
    print("="*70)
    print("\nFeatures:")
    print("  ✓ Interactive folder and file selection")
    print("  ✓ Excel report parsing with validation")
    print("  ✓ Dynamic data filtering by row selection")
    print("  ✓ Multiple chart types (Pie, Bar, Heatmap, Scatter)")
    print("  ✓ High-quality visualizations with Plotly")
    print("  ✓ Responsive web interface")
    print("="*70)
    print("\nPress CTRL+C to quit\n")
    
    # Create and run application
    try:
        app = PathogenReportAnalyzerApp()
        app.run(host=args.host, port=args.port, debug=args.debug)
    except KeyboardInterrupt:
        print("\n\n" + "="*70)
        print("Server stopped by user")
        print("="*70)
        sys.exit(0)
    except Exception as e:
        print(f"\n❌ Error starting server: {e}")
        traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()