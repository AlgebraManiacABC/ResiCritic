def flatten_nested_list(nested: list[list[str]]):
    return [item for sublist in nested for item in sublist]

def unique(lst: list):
    return list(set(lst))

def unique_flattened_nested_str_list(nested: list[list[str]]):
    return [val for val in set(flatten_nested_list(nested)) if val != ""]

def remove_empty(lst: list):
    return [val for val in lst if val]
