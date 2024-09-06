import pandas as pd

def format_excel(writer, df, sheet_name, header_format=None, cell_format=None):
    workbook = writer.book
    worksheet = writer.sheets[sheet_name]

    if header_format:
        header_fmt = workbook.add_format(header_format)
        for col_num, value in enumerate(df.columns.values):
            worksheet.write(0, col_num, value, header_fmt)

    if cell_format:
        cell_fmt = workbook.add_format(cell_format)
        for row in range(1, len(df) + 1):
            for col in range(len(df.columns)):
                worksheet.write(row, col, df.iloc[row - 1, col], cell_fmt)

def csv_to_excel(csv_file_path, excel_file_path, sheet_name='Sheet1'):
    df = pd.read_csv(csv_file_path)

    with pd.ExcelWriter(excel_file_path, engine='xlsxwriter', options={'nan_inf_to_errors': True}) as writer:

        header_format = {
            'bold': True,
            'text_wrap': True,
            'valign': 'top',
            'fg_color': '#D7E4BC',
            'border': 1
        }
        cell_format = {
            'border': 1
        }
        format_excel(writer, df, sheet_name, header_format, cell_format)


if __name__ == '__main__':

    csv_file = '/Users/josediogomoura/Desktop/BioFago/github/data/output/phages/phages_info.csv'
    excel_file = '/Users/josediogomoura/Desktop/BioFago/github/data/output/phages/output_excel_file.xlsx'

    csv_to_excel(csv_file, excel_file)
