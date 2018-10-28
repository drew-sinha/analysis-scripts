import datetime
import collections

from IPython import display
import pandas

# Take a list and unfold it all the way
def flatten_list(my_list, to_level=-1, this_level=0):
    flat_list = []

    for my_item in my_list:
        if (not isinstance(my_item, collections.Iterable)) or (type(my_item) is str) or ((to_level is not -1) and (this_level is to_level)):
            flat_list.append(my_item)
        else:
            flat_list.extend(flatten_list(my_item,to_level,this_level+1))
    return flat_list

def extract_datetime_fromstr(time_str):
    return datetime.datetime.strptime(time_str,'%Y-%m-%dt%H%M')

def unique_items(seq):  # Credit to Peterbe
    seen = set()
    return [x for x in seq if x not in seen and not seen.add(x)]

def in_range(arr,low,high):
    return (arr >= low) & (arr <= high) & (~np.isnan(arr))


def print_table(array_table, column_names=None, row_names=None, optional_formatting=None):
    if optional_formatting is not None:
        array_table = [[optional_formatting.format(val) for val in row] for row in array_table]

    # if row_names:
    #     new_table = []
    #     for row_name, array_row in zip(row_names, array_table):
    #         base = [row_name]
    #         base.extend(array_row)
    #         new_table.append(base)
    #     array_table = new_table
    #
    # if column_names:
    #     column_names = [column_names]
    #     if row_name:
    #         temp = [None]
    #         temp.extend(column_names)
    #         column_names = temp
    #
    #     column_names.append(array_table)
    #     array_table = column_names

    formatting_df = pandas.DataFrame(data=array_table, columns=column_names, index=row_names)
    display.display(formatting_df)
