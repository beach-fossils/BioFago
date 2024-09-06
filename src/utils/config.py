import yaml
from pathlib import Path


class Config:
    def __init__(self, config_file: str = 'config.yaml'):
        self.config_file = Path(config_file)
        self.load_config()

    def load_config(self):
        try:
            if not self.config_file.exists():
                raise FileNotFoundError(f"Config file not found: {self.config_file}")

            with open(self.config_file, 'r') as file:
                config = yaml.safe_load(file)

            self.genomes_folder = config.get('genomes_folder', '')
            if not self.genomes_folder:
                raise ValueError("genomes_folder not specified in config file")

            self.keep_sequence_loci = config.get('keep_sequence_loci', False)
            self.threshold_species = config.get('threshold_species', 0.95)
            self.output_folder = config.get('output_folder', '')
            self.log_level = config.get('log_level', 'INFO')


        except yaml.YAMLError as e:
            print(f"Error parsing YAML file: {e}")

        except Exception as e:
            print(f"Error loading configuration: {e}")

    def __str__(self):
        return f"Config(genomes_folder='{self.genomes_folder}', keep_sequence_loci={self.keep_sequence_loci}, threshold_species={self.threshold_species})"

# Usage example:
# config = Config()
# print(config.genomes_folder)
