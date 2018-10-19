import datetime
import collections

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
