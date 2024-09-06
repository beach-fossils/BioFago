import logging
from pathlib import Path


def setup_logging(log_file: Path, log_level: str = 'INFO'):
    """
    Set up logging configuration.

    :param log_file: Path to the log file
    :param log_level: Logging level (default: INFO)
    """
    numeric_level = getattr(logging, log_level.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError(f'Invalid log level: {log_level}')

    logging.basicConfig(
        level=numeric_level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
        filename=str(log_file),
        filemode='w'
    )

    # Add console handler to display logs in console as well
    console = logging.StreamHandler()
    console.setLevel(numeric_level)
    formatter = logging.Formatter('%(name)-12s: %(levelname)-8s %(message)s')
    console.setFormatter(formatter)
    logging.getLogger('').addHandler(console)

    logging.info(f"Logging initialized. Log file: {log_file}")

# Usage example:
# setup_logging(Path('path/to/your/logfile.log'), 'INFO')