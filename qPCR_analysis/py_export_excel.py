import pandas as pd
import numpy as np
import xlsxwriter
import itertools

if df.shape != var.shape:
    raise ValueError('values and error bars tables ar not matching in size')

# the first element is empty because when you have a df eg. with 1 column (df.shape[1] == 1)  you want the letter A and not B
col_letter = ['', 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z', 'AA', 'AB', 'AC', 'AD', 'AE', 'AF', 'AG', 'AH', 'AI', 'AJ', 'AK', 'AL', 'AM', 'AN', 'AO', 'AP', 'AQ', 'AR', 'AS', 'AT', 'AU', 'AV', 'AW', 'AX', 'AY', 'AZ']

def col_to_merge(init_row, init_col, length):
    """returns a range of cells like A2:F2"""
    if init_col < 1:
        raise ValueError("Index < 1 not allowed for columns, the function assumes index 1 for column A")
    return col_letter[init_col] + str(init_row) + ':' + col_letter[init_col-1+length] + str(init_row)


# # Paste the df on excel


def place_in_excel(init_row, init_col, df, title = ''):
    if title != '': # if the function was called with a title, than plot it as a merged cell
        worksheet.merge_range(col_to_merge(init_row, init_col, len(df.columns)), title )
   

    #Add the column names
    for i in range(0,len(df.columns)):
                        # row, col, argument
        worksheet.write(init_row, init_col+i-1, df.columns[i])
        

    #Add the content
    for index, row in df.iterrows():
        for i in range(0,len(row)):
            worksheet.write(init_row+1+index, init_col+i-1 ,row[i])

# # Multiple plot Function

def multiple_plot(df):
    """generates multiple plot, one for each gene of interest"""
    for i in range(0, df.shape[0]): #number of rows
        chart = workbook.add_chart({'type' : 'column'})
        
        chart.add_series({'name' : '=Sheet1!A{}'.format(str(i+3)), # the name of the genes are in the column A, starting from row 3
                    'values': '=Sheet1!B{0}:{1}{0}'.format(
                                                        str(i+3),
                                                        col_letter[df.shape[1]] # define ending columns based on the number of columns of the df
                          ),
                    'categories' : '=Sheet1!B2:{}2'.format(col_letter[df.shape[1]]),
                    'y_error_bars': {
                        'type':         'custom',
                        'plus_values':  '=Sheet1!B{0}:{1}{0}'.format(
                                                                    str(i+2+df.shape[0]+3), #df.shape[0]+2 is the offset for placing the errors table
                                                                    col_letter[df.shape[1]]
                        ),
                        'minus_values': '=Sheet1!B{0}:{1}{0}'.format(
                                                                    str(i+2+df.shape[0]+3), #df.shape[0]+2 is the offset for placing the errors table
                                                                    col_letter[df.shape[1]]
                        )}
        })
        
        # position each chart at a distance of 2 columns from the end of the df and 20 rows from each others 
        worksheet.insert_chart(col_letter[df.shape[1]+3] + str(i*20+1), chart)
        
        # print(df.iloc[i]) #implicit indexing based on the actual structure and order of the df


# # Generate the excel file 

def generate_xlsx(value_df, variance_df):

    place_in_excel(1, 1, value_df, title = "dCT values")

    place_in_excel(value_df.shape[0]+3, 1, variance_df, title = "Std error")

    multiple_plot(value_df)

workbook = xlsxwriter.Workbook('result.xlsx')
worksheet = workbook.add_worksheet()

generate_xlsx(df, var)

workbook.close()  



