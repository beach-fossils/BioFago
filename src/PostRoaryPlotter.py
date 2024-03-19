import json
from matplotlib.colors import ListedColormap
import re
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import logging


def setup_logging(self, log_file):
    logging.basicConfig(
        filename=log_file,
        filemode='a',
        format='%(asctime)s - %(levelname)s - %(message)s',
        level=logging.INFO
    )


logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')


class GenePresenceAbsencePlotter:
    def __init__(self, path_csv):
        self.path_csv = path_csv
        self.df_binary = None
        self.df_selected = None
        self.pattern_counts = None
        self.organism_names = {}
        self.df_final = None  # This will hold the final dataframe

    def process_strain_names(self):
        new_columns = []
        for col in self.df_binary.columns:
            # This regex matches the accession numbers including their version numbers
            # It correctly captures the prefix (GCA or GCF), the numeric part, and the version number after the dot
            match = re.match(r'^(GCA|GCF)\.(\d+)\.(\d+)', col)
            if match:
                # Format the simplified name to include the underscore and retain the version number as specified
                simplified_name = f"{match.group(1)}_{match.group(2)}.{match.group(3)}"
                new_columns.append(simplified_name)
            else:
                new_columns.append(col)  # Keep the original column name if no match is found
        self.df_binary.columns = new_columns

    def load_and_process_data(self):
        df = pd.read_csv(self.path_csv, sep=',', index_col=0)
        index_avg_group_size = df.columns.get_loc('Avg group size nuc')
        self.df_binary = df.iloc[:, index_avg_group_size + 1:].notna().astype(int)
        self.process_strain_names()

    def identify_unique_patterns(self):
        # Transpose df_binary so that genes are columns and strains are rows
        df_transposed = self.df_binary.transpose()

        # Create a summary representation for each strain
        pattern_strings = df_transposed.apply(lambda x: ''.join(x.astype(str)), axis=1)
        self.pattern_counts = pattern_strings.value_counts()

        # Select representative strains
        representative_strains = pattern_strings.drop_duplicates()
        self.df_selected = df_transposed.loc[representative_strains.index]

    def plot_heatmap(self):
        # Create a custom color map for the heatmap
        cmap = ListedColormap(['#ffffd9', '#41b6c4'])  # Yellow for absence (0), Blue for presence (1)

        # Plotting the heatmap
        plt.figure(figsize=(20, 10))
        ax = sns.heatmap(self.df_selected, cmap=cmap, cbar=True, xticklabels=True, yticklabels=True,
                         cbar_kws={'ticks': [0, 1], 'label': 'Gene Presence/Absence'})

        # Customize the color bar to show "Absent" and "Present" labels
        cbar = ax.collections[0].colorbar
        cbar.set_ticklabels(['Absent', 'Present'])

        # Adjusting the annotation process to match strain names with their counts
        for i, (index_value, row) in enumerate(self.df_selected.iterrows()):
            pattern = ''.join(row.astype(str))
            count = self.pattern_counts[pattern]
            plt.text(0.5, i + 0.5, f'Count: {count}', ha='center', va='center', color='black', fontsize=7)

        plt.title('Heatmap of Gene Presence/Absence Across Selected Strains')
        plt.xlabel('Genes')
        plt.ylabel('Selected Strains (Unique Patterns)')
        plt.xticks(rotation=90)  # Rotate gene labels for better readability
        plt.show()

    def get_species_names(self, ncbi_jsonl_path):
        with open(ncbi_jsonl_path, 'r') as file:
            for line in file:
                record = json.loads(line)
                # Accession is directly available
                accession = record.get('accession', 'Unknown')
                # Organism name is nested within several dictionaries
                organism_name = record.get('assemblyInfo', {}).get('biosample', {}).get('description', {}).get(
                    'organism', {}).get('organismName', 'Unknown')
                self.organism_names[accession] = organism_name

        return self.organism_names

    def generate_final_csv(self, output_csv_path):
        # First, ensure the organism_names dictionary is populated
        if not self.organism_names:
            return

        # Prepare the data for CSV output
        output_data = {
            "accession": [],
            "organismName": []
        }

        # Add gene columns to the output data
        for gene in self.df_binary.index:
            output_data[gene] = []

        # Iterate over each column (strain) in the binary matrix
        for col in self.df_binary.columns:
            # Add the accession number and organism name
            output_data["accession"].append(col)
            output_data["organismName"].append(self.organism_names.get(col, "Unknown"))

            # Add the presence/absence data for each gene
            for gene in self.df_binary.index:
                output_data[gene].append(self.df_binary.loc[gene, col])

        # Convert the output data to a DataFrame and save as CSV
        output_df = pd.DataFrame(output_data)
        output_df.to_csv(output_csv_path, index=False)

        logging.info(f"CSV file saved to {output_csv_path}")

        return output_df

    def make_final_csv(self, output_csv_path, json_path):
        try:
            self.get_species_names(
                json_path)

            self.generate_final_csv(
                output_csv_path)

        except Exception as e:
            logging.error(f"An error occurred: {e}")


class GenePresenceAbsenceCluster:
    def __init__(self, data_path):
        self.data_path = data_path
        self.data = None
        self.grouped_data = None

    def load_data(self):
        self.data = pd.read_csv(self.data_path)
        self.data.drop('accession', axis=1, inplace=True)

    def process_data(self):
        self.data['pattern'] = self.data.drop('organismName', axis=1).apply(
            lambda row: ','.join(row.values.astype(str)), axis=1)
        self.grouped_data = self.data.groupby(['organismName', 'pattern']).first().reset_index()
        self.grouped_data.drop('pattern', axis=1, inplace=True)
        self.grouped_data.set_index('organismName', inplace=True)

    def plot_heatmap(self):
        sns.clustermap(self.grouped_data, cmap='viridis', linewidths=.5, annot=False,
                       cbar_kws={'label': 'Gene Presence/Absence'})
        plt.title('Clustermap of Gene Presence/Absence Across Species')
        plt.xlabel('Genes')
        plt.ylabel('Species')
        plt.show()

    def process_and_plot(self):
        try:
            self.load_data()
            self.process_data()
            self.plot_heatmap()
        except Exception as e:
            logging.error(f"An error occurred: {e}")


class GenePresenceAbsencePivot:
    def __init__(self, data_path):
        self.data_path = data_path
        self.data = None
        self.pattern_counts = None
        self.pivot_table = None

    def load_data(self):
        self.data = pd.read_csv(self.data_path)
        self.data.drop('accession', axis=1, inplace=True)

    def process_data(self):
        self.data['pattern'] = self.data.drop('organismName', axis=1).apply(
            lambda row: ','.join(row.values.astype(str)), axis=1)
        self.pattern_counts = self.data.groupby(['organismName', 'pattern']).size().reset_index(name='count')

    def create_pivot_table(self):
        self.pivot_table = self.pattern_counts.pivot(index='organismName', columns='pattern', values='count').fillna(0)

    def plot_heatmap(self):
        plt.figure(figsize=(20, 10))
        sns.heatmap(self.pivot_table, annot=True, fmt="g", cmap='mako', cbar_kws={'label': 'Count of Unique Patterns'})
        plt.title('Heatmap of Unique Gene Presence/Absence Patterns Across Species')
        plt.xlabel('Unique Patterns of Gene Presence/Absence')
        plt.ylabel('Species')
        plt.xticks(rotation=45, ha='right')  # Rotate the x labels for better readability
        plt.show()

    def process_and_plot(self):
        self.load_data()
        self.process_data()
        self.create_pivot_table()
        self.plot_heatmap()


