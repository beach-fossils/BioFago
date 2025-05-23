import logging
import re
import csv
import os
from typing import Dict, List, Tuple, Set, Union
from statistics import mean, stdev

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
matrix_path = os.path.join(SCRIPT_DIR, '..', '..', 'reference_crispr', 'clades_groups_correlation.csv')

class CRISPRCladeClassifier:
    def __init__(self, matrix_path: str):
        self.clade_matrix = self.load_clade_matrix(matrix_path)
        self.clade_data = {
            "Widely Prevalent (WP)": {
                "Ia": ((28, 36), (31, 35), 5),
                "Ib": ((27, 36), (21, 26), 5),
                "II": ((12, 19), (27, 35), 5)
            },
            "Western North American (WNA)": {"III": ((44, 110), (25, 49), 5)},
            "Eastern North American (ENA)": {"IV": ((53, 95), (25, 58), 5)},
            "B-Group": {
                "I": ((29, 31), (22, 23), 5),
                "II": ((18, 22), (49, 49), 5),
                "III": ((18, 23), (36, 39), 5),
                "IV": ((42, 42), (32, 32), 5),
                "ATCC-BAA2158": ((53, 53), (38, 38), 5)
            },
            "Rubus": {
                "I": ((20, 55), (16, 40), (3, 6)),
                "II": ((20, 41), (16, 38), (3, 6)),
                "III": ((48, 55), (31, 40), (3, 6))
            }
        }

    @staticmethod
    def load_clade_matrix(matrix_path: str) -> List[Dict[str, str]]:
        with open(matrix_path, 'r') as file:
            return list(csv.DictReader(file))

    def parse_crispr_genotype(self, genotype_str: str) -> Dict[str, Dict[str, Union[str, float]]]:
        logging.info(f"Parsing genotype string: {genotype_str}")
        crr_data = {}

        if not genotype_str or genotype_str.lower() == 'none':
            logging.warning("Empty or None genotype string")
            return crr_data

        # Remove "CRR: " prefix if present
        genotype_str = re.sub(r'^CRR:\s*', '', genotype_str)

        # Split by semicolon and process each part
        for crr_info in genotype_str.split(';'):
            crr_info = crr_info.strip()
            # Match both formats
            pattern = r'CRR(\d+)\s*-?\s*(?:Group:\s*)?(.+?),\s*Score:\s*([\d.]+)'
            match = re.match(pattern, crr_info)

            if match:
                crr_num, info, score = match.groups()
                # Add "Group: " prefix if not present
                if not info.startswith('Group:'):
                    info = f"Group: {info}"
                crr_data[f'CRR{crr_num}'] = {'info': info, 'score': float(score)}
                logging.info(f"Parsed CRR{crr_num}: info='{info}', score={score}")
            else:
                logging.warning(f"Could not parse CRR info: {crr_info}")

        return crr_data

    def extract_groups(self, crr_info: str) -> Tuple[str, str]:
        # WP cases
        wp_match = re.search(r'WP,\s*Subgroup:\s*Group\s+(Ia \+ Ib|Ia \+ II|II|Ib)', crr_info)
        if wp_match:
            return "WP", wp_match.group(1)

        # Other specific groups
        group_match = re.search(r'Group (IH 3-1|E\.pyrifoliae DSM 12163)', crr_info)
        if group_match:
            return group_match.group(1), ""

        # ENA and WNA
        if "Group IV" in crr_info:
            return "ENA", "IV"
        if "Group III" in crr_info:
            return "WNA", "III"

        # Rubus cases
        rubus_match = re.search(r'Rubus,\s*Subgroup:\s*Rubus\s+(I|II|III)', crr_info)
        if rubus_match:
            return "Rubus", rubus_match.group(1)

        # B-group cases
        b_group_match = re.search(r'B-group(?:[,-]\s*(I{1,4}|ATCC(?:-BAA2158)?))?', crr_info, re.IGNORECASE)
        if b_group_match:
            return "B-group", b_group_match.group(1) if b_group_match.group(1) else ""

        return "", ""

    def determine_clade(self, genotype_str: str, spacer_counts: Dict[str, int]) -> Tuple[
        str, float, float, float, str, str]:
        def parse_crr_info(crr_info: str) -> Tuple[str, str, str]:
            # Modified to better catch WP groups
            group_match = re.search(r'Group:\s*([\w.-]+(?:\s+[\w.-]+)*(?:\s*\([^)]*\))*)', crr_info)
            subgroup_match = re.search(r'Subgroup:\s*(?:Group\s+)?([^,>]+)', crr_info)  # Made 'Group' optional
            specific_match = re.search(r'>\s*([^,]+)', crr_info)

            group = group_match.group(1) if group_match else ""
            subgroup = subgroup_match.group(1) if subgroup_match else ""
            specific = specific_match.group(1) if specific_match else ""
            return group, subgroup, specific

        def extract_wp_subgroups(subgroup: str) -> List[str]:
            # More flexible pattern to catch variations
            return re.findall(r'(Ia|Ib|II)(?:\s*\+\s*)?', subgroup)

        def extract_subgroup(subgroup: str) -> str:
            subgroup_match = re.search(r'(Rubus\s+[IV]+|B-group\s+[IV]+|Group\s+[IV]+|[IV]+)', subgroup)
            return subgroup_match.group(1) if subgroup_match else ""

        logging.info(f"Determining clade for genotype: {genotype_str}")
        logging.info(f"With spacer counts: {spacer_counts}")

        if not genotype_str or genotype_str.lower() == 'none':
            return "Unknown", 0.0, 0.0, 0.0, "Very Low", ""

        # Remove "CRR: " prefix if present
        genotype_str = re.sub(r'^CRR:\s*', '', genotype_str)
        crr_data = self.parse_crispr_genotype(genotype_str)
        logging.info(f"Parsed CRR data: {crr_data}")

        crr1_info = crr_data.get('CRR1', {}).get('info', '')
        crr2_info = crr_data.get('CRR2', {}).get('info', '')
        crr1_score = crr_data.get('CRR1', {}).get('score', 0)
        crr2_score = crr_data.get('CRR2', {}).get('score', 0)

        logging.info(f"CRR1 info: {crr1_info}, score: {crr1_score}")
        logging.info(f"CRR2 info: {crr2_info}, score: {crr2_score}")

        crr1_group, crr1_subgroup, crr1_specific = parse_crr_info(crr1_info)
        crr2_group, crr2_subgroup, crr2_specific = parse_crr_info(crr2_info)

        # WP cases first
        if "WP" in crr1_group and "WP" in crr2_group:
            clade = "Widely Prevalent (WP)"
            crr1_subgroups = extract_wp_subgroups(crr1_subgroup)
            crr2_subgroups = extract_wp_subgroups(crr2_subgroup)

            common_subgroups = set(crr1_subgroups) & set(crr2_subgroups)
            if common_subgroups:
                subgroup = "+".join(sorted(common_subgroups))
            else:
                subgroup = f"{'+'.join(sorted(crr1_subgroups))}/{'+'.join(sorted(crr2_subgroups))}"

        # ENA and WNA cases
        elif all(("ENA" in g or "Group IV" in g) for g in [crr1_group, crr2_group]):
            clade = "Eastern North American (ENA)"
            subgroup = "IV"
        elif all(("WNA" in g or "Group III" in g) for g in [crr1_group, crr2_group]):
            clade = "Western North American (WNA)"
            subgroup = "III"

        # Rubus cases
        elif all("Rubus" in g for g in [crr1_group, crr2_group]):
            clade = "Rubus"
            rubus_subgroup = extract_subgroup(crr1_subgroup) or extract_subgroup(crr2_subgroup)
            subgroup = rubus_subgroup.replace("Rubus ", "") if rubus_subgroup else "I"

        # B-group cases
        elif all("B-group" in g for g in [crr1_group, crr2_group]):
            clade = "B-Group"
            b_group_subgroup = extract_subgroup(crr1_subgroup) or extract_subgroup(crr2_subgroup)
            subgroup = b_group_subgroup.replace("B-group ", "") if b_group_subgroup else "II"

        # Special cases
        elif all("E.pyrifoliae DSM 12163" in g for g in [crr1_group, crr2_group]):
            clade = "E.pyrifoliae DSM 12163"
            subgroup = ""
        elif all("IH 3-1" in g for g in [crr1_group, crr2_group]):
            clade = "IH 3-1"
            subgroup = ""
        else:
            return "Unknown", 0.0, 0.0, 0.0, "Very Low", ""

        # Calculate scores
        genotype_score = (crr1_score + crr2_score) / 2

        # Calculate spacer score if we have matching clade data
        if clade in self.clade_data and subgroup in self.clade_data[clade]:
            crr1_range, crr2_range, crr4_range = self.clade_data[clade][subgroup]
            spacer_score = self.calculate_spacer_score(
                spacer_counts.get('CRR1', 0),
                spacer_counts.get('CRR2', 0),
                spacer_counts.get('CRR4', 0),
                crr1_range, crr2_range, crr4_range
            )
        else:
            spacer_score = 0

        final_confidence_score = 0.7 * genotype_score + 0.3 * spacer_score
        confidence_level = self.determine_confidence_level(final_confidence_score)

        return clade, final_confidence_score, spacer_score, genotype_score, confidence_level, subgroup

    @staticmethod
    def determine_confidence_level(score: float) -> str:
        if score >= 95:
            return "Very High"
        elif score >= 85:
            return "High"
        elif score >= 70:
            return "Moderate"
        elif score >= 50:
            return "Low"
        else:
            return "Very Low"

    @staticmethod
    def calculate_spacer_score(crr1_spacers: int, crr2_spacers: int, crr4_spacers: int,
                               crr1_range: Tuple[int, int], crr2_range: Tuple[int, int],
                               crr4_range: Union[int, Tuple[int, int]]) -> float:
        if crr1_spacers == 0 and crr2_spacers == 0 and crr4_spacers == 0:
            return 0

        scores = []
        # Score CRR1
        if crr1_range[0] <= crr1_spacers <= crr1_range[1]:
            scores.append(100)
        else:
            diff = min(abs(crr1_spacers - crr1_range[0]), abs(crr1_spacers - crr1_range[1]))
            scores.append(max(0, 100 - (diff * 5)))

        # Score CRR2
        if crr2_range[0] <= crr2_spacers <= crr2_range[1]:
            scores.append(100)
        else:
            diff = min(abs(crr2_spacers - crr2_range[0]), abs(crr2_spacers - crr2_range[1]))
            scores.append(max(0, 100 - (diff * 5)))

        # Score CRR4
        if isinstance(crr4_range, tuple):
            if crr4_range[0] <= crr4_spacers <= crr4_range[1]:
                scores.append(100)
            else:
                diff = min(abs(crr4_spacers - crr4_range[0]), abs(crr4_spacers - crr4_range[1]))
                scores.append(max(0, 100 - (diff * 5)))
        else:
            scores.append(100 if crr4_spacers == crr4_range else max(0, 100 - abs(crr4_spacers - crr4_range) * 5))

        return sum(scores) / len(scores)

    @staticmethod
    def calculate_combined_score(spacer_score: float, genotype_score: float) -> float:
        return (spacer_score * 0.3 + genotype_score * 0.7)  # Give more weight to genotype score

    @staticmethod
    def calculate_final_confidence_score(spacer_score: float, genotype_score: float) -> float:
        return (spacer_score * 0.2 + genotype_score * 0.8)  # Give even more weight to genotype score


    def determine_wp_subgroup(self, crr1_spacers: int, crr2_spacers: int, crr4_spacers: int) -> Tuple[str, float]:
        wp_ranges = self.clade_data["Widely Prevalent (WP)"]

        scores = {}
        for subgroup, (crr1_range, crr2_range, crr4_count) in wp_ranges.items():
            crr1_score = self.calculate_range_score(crr1_spacers, range(crr1_range[0], crr1_range[1] + 1))
            crr2_score = self.calculate_range_score(crr2_spacers, range(crr2_range[0], crr2_range[1] + 1))
            crr4_score = self.calculate_range_score(crr4_spacers, [crr4_count])
            score = (crr1_score + crr2_score + crr4_score) / 3
            scores[subgroup] = score

            print(f"Debug - WP subgroup {subgroup} scores: CRR1: {crr1_score}, CRR2: {crr2_score}, CRR4: {crr4_score}, Mean: {score}")

        best_subgroup = max(scores, key=scores.get)
        return best_subgroup, scores[best_subgroup]

    @staticmethod
    def determine_confidence_level(score: float) -> str:
        if score >= 99:
            return "Excellent"
        elif score >= 95:
            return "Very High"
        elif score >= 85:
            return "High"
        elif score >= 70:
            return "Moderate"
        elif score >= 60:
            return "Low"
        else:
            return "Very Low"

    @staticmethod
    def calculate_range_score(value: int, range_obj: Union[range, List[int]]) -> float:
        if isinstance(range_obj, range):
            if value in range_obj:
                return 100.0

            min_val, max_val = range_obj.start, range_obj.stop - 1
            if value < min_val:
                distance = min_val - value
                return max(0, 100 - (distance * 10))
            else:
                distance = value - max_val
                return max(0, 100 - (distance * 5))
        else:
            if value == range_obj[0]:
                return 100.0
            distance = abs(value - range_obj[0])
            return max(0, 100 - (distance * 10))

    @staticmethod
    def parse_range(range_str: str) -> Union[range, List[int]]:
        if '-' in range_str:
            start, end = map(int, range_str.split('-'))
            return range(start, end + 1)
        else:
            return [int(range_str)]


    @staticmethod
    def calculate_genotype_score(crr1_score: float, crr2_score: float, common_groups: Set[str]) -> Tuple[float, float]:
        if crr1_score or crr2_score:
            if common_groups:
                genotype_score = (crr1_score + crr2_score) / 2
                genotype_uncertainty = abs(crr1_score - crr2_score) / 2
            else:
                genotype_score = min(crr1_score, crr2_score) * 0.7
                genotype_uncertainty = abs(crr1_score - crr2_score) / 2
        else:
            genotype_score = 0
            genotype_uncertainty = 0
        return genotype_score, genotype_uncertainty


def parse_spacer_counts(spacer_str: str) -> Dict[str, int]:
    """Parse CRISPR spacer counts from log output."""
    logging.info(f"Parsing spacer counts from: {spacer_str}")

    if not spacer_str or spacer_str == "No results":
        logging.warning("No spacer string provided or 'No results' found")
        return {}

    counts = {}

    # First try to find CRRFinder's detailed counts
    if isinstance(spacer_str, str):  # Add type check
        lines = spacer_str.split('\n')
        for line in lines:
            match = re.search(r'Found (\d+) (CRR[124])', line)
            if match:
                count, crr_type = match.groups()
                counts[crr_type] = int(count)

        # If no matches found, try CRISPRAnalyzer format
        if not counts and ': ' in spacer_str:
            parts = [p.strip() for p in spacer_str.split(',')]
            for part in parts:
                if ': ' in part and 'Total:' not in part:
                    crr, count = part.split(': ')
                    try:
                        counts[crr] = int(count)
                    except ValueError:
                        logging.error(f"Could not convert count to integer in part: {part}")
                        continue

    logging.info(f"Parsed spacer counts: {counts}")
    return counts

def process_test_case(classifier: CRISPRCladeClassifier, case: Union[Tuple[str, str], str]) -> Dict[str, Union[str, float]]:
    if isinstance(case, tuple) and len(case) == 2:
        genotype, spacers = case
    else:
        genotype = case
        spacers = ""

    if (not genotype or genotype.lower() == 'none') and not spacers:
        return {"clade": "Insufficient data to determine clade."}

    spacer_counts = parse_spacer_counts(spacers)
    clade, score, spacer_score, genotype_score, confidence_level, subgroup = classifier.determine_clade(genotype, spacer_counts)

    return {
        "clade": f"{clade} {subgroup}".strip() if clade != "Unknown" else clade,
        "confidence_score": score,
        "confidence_level": confidence_level,
        "spacer_match_score": spacer_score,
        "genotype_score": genotype_score,
        "input_genotype": genotype,
        "input_spacers": spacers
    }

def run_classifier(matrix_path: str, test_cases: List[Union[Tuple[str, str], str]]) -> List[Tuple[int, Dict[str, Union[str, float]]]]:
    classifier = CRISPRCladeClassifier(matrix_path)
    return [(i, process_test_case(classifier, case)) for i, case in enumerate(test_cases, 1)]


def process_csv(input_file: str, output_file: str, classifier: CRISPRCladeClassifier):
    with open(input_file, 'r') as infile, open(output_file, 'w', newline='') as outfile:
        reader = csv.DictReader(infile)
        fieldnames = reader.fieldnames + ['Clade', 'Confidence_Score', 'Confidence_Level']
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        writer.writeheader()

        for row in reader:
            genotype = row['crispr_genotype_extra_info']
            spacers = row['crispr_spacers_count']

            spacer_counts = parse_spacer_counts(spacers)
            clade, score, spacer_score, genotype_score, confidence_level, subgroup = classifier.determine_clade(
                genotype, spacer_counts)

            row['Clade'] = f"{clade} {subgroup}".strip() if clade != "Unknown" else clade
            row['Confidence_Score'] = f"{score:.2f}"
            row['Confidence_Level'] = confidence_level

            writer.writerow(row)



def main():
    test_cases = [
        (
        "CRR2 - Group: WP, Subgroup: Group Ia + II > a_Ea3a-H1-II, Score: 100.00; CRR1 - Group: WP, Subgroup: Group Ia + Ib > E_3477-1, Score: 92.72",
        "CRR1: 27, CRR2: 34, CRR4: 5, Total: 66"),
        (
        "CRR: CRR4 - Group: Group Spiraeoideae, Subgroup: alpha, Score: 100.00; CRR2 - Group: WP, Subgroup: Group Ia + II > a_Ea3a-H1-II, Score: 100.00; CRR1 - Group: WP, Subgroup: Group Ia + Ib > A, Score: 100.00",
        "CRR1: 27, CRR2: 34, CRR4: 5, Total: 66"),
        ("CRR2 - Group: WP, Subgroup: Group II, Score: 100.00; CRR1 - B-group, Subgroup: B-group II, Score: 95.00",
         "CRR1: 20, CRR2: 49, CRR4: 5, Total: 74"),
        ("CRR2 - Group: WP, Subgroup: Group Ia + II, Score: 100.00; CRR1 - Group: WP, Subgroup: Group II, Score: 92.72",
         "CRR1: 15, CRR2: 30, CRR4: 5, Total: 50"),
        ("CRR2 - Group: WP, Subgroup: Group Ia, Score: 100.00; CRR1 - Group: WP, Subgroup: Group Ib, Score: 92.72",
         "CRR1: 30, CRR2: 25, CRR4: 5, Total: 60"),
        ("CRR2 - Group: WNA (Group III), Score: 100.00; CRR1 - Group: ENA (Group IV), Score: 92.72",
         "CRR1: 60, CRR2: 30, CRR4: 5, Total: 95"),
        (
        "CRR: CRR4 - Group: Group Spiraeoideae, Subgroup: alpha, Score: 100.00; CRR2 - Group: WP, Subgroup: Group Ia + II > a_57669, Score: 100.00; CRR1 - Group: WP, Subgroup: Group Ia + Ib > D_57669_delete, Score: 100.00",
        "CRR1: 27, CRR2: 34, CRR4: 5, Total: 66"),
        (
        "CRR4 - Group: Group Spiraeoideae, Subgroup: alpha, Score: 100.00; CRR2 - Group: WP, Subgroup: Group Ia+Ib, Score: 97.06; CRR1 - Group: WP, Subgroup: Group Ia+Ib, Score: 88.24",
        "CRR1: 27, CRR2: 34, CRR4: 5, Total: 66"),
        (
        "CRR: CRR4 - Group: Group Rubus, Subgroup: Ea1-97, Score: 100.00; CRR2 - Group: Rubus, Subgroup: Rubus I > Ea1-97, Score: 99.00; CRR1 - Group: Rubus, Subgroup: Rubus I > Ea1-97, Score: 96.00",
        "CRR1: 60, CRR2: 30, CRR4: 6, Total: 96"),
        (
        "CRR: CRR4 - Group: Group Spiraeoideae, Subgroup: Ea01-03, Score: 100.00; CRR2 - Group: B-group, Subgroup: B-group IV > Ea01-03, Score: 100.00; CRR1 - Group: B-group, Subgroup: B-group IV > Ea01-03, Score: 100.00 ",
        ""),
        ("None", "CRR1: 60, CRR2: 55, CRR4: 5, Total: 95"),
        (
        "CRR: CRR4 - Group: Group Spiraeoideae, Subgroup: alpha, Score: 100.00; CRR2 - Group: ENA (Group IV), Subgroup: MLI200-18, Score: 100.00; CRR1 - Group: ENA (Group IV), Subgroup: MLI200-18, Score: 100.00",
        ""),
        (
        "CRR: CRR4 - Group: Group Spiraeoideae, Subgroup: alpha, Score: 100.00; CRR2 - Group: WNA (Group III), Subgroup: 1280, Score: 100.00; CRR1 - Group: WNA (Group III), Subgroup: 7-3, Score: 100.00",
        ""),
        (
        "CRR: CRR4 - Group: Group Spiraeoideae, Subgroup: alpha, Score: 100.00; CRR2 - Group: E.pyrifoliae DSM 12163, Subgroup: E.pyrifoliae DSM 12163, Score: 100.00; CRR1 - Group: E.pyrifoliae DSM 12163, Subgroup: E.pyrifoliae DSM 12163, Score: 100.00",
        ""),
        (
        "CRR: CRR4 - Group: Group Spiraeoideae, Subgroup: alpha, Score: 100.00; CRR2 - Group: IH 3-1, Subgroup: IH 3-1, Score: 100.00; CRR1 - Group: IH 3-1, Subgroup: IH 3-1, Score: 100.00",
        ""),
    ]

    results = run_classifier(matrix_path, test_cases)

    for i, result in results:
        print(f"Test case {i} clade: {result['clade']}")
        print(f"Confidence score: {result['confidence_score']:.2f}")
        print(f"Confidence level: {result['confidence_level']}")
        print(f"Spacer match score: {result['spacer_match_score']:.2f}")
        print(f"Genotype score: {result['genotype_score']:.2f}")
        print(f"Input: {result['input_genotype']}")
        print(f"Spacers: {result['input_spacers']}\n")


if __name__ == '__main__':
    main()