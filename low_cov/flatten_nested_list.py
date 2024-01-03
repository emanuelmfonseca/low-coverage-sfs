import itertools

def flatten_nested_list(nested_list, operation):
    """
    Flatten a nested list by concatenating or multiplying elements.
    
    Args:
        nested_list (list of lists): The nested list to be flattened.
        operation (str): The operation to apply between elements ('+' for concatenation or '*' for multiplication).
    
    Returns:
        list: A flattened list resulting from applying the operation to elements.
    """
    if len(nested_list) == 1:
        # If there's only one sublist, return it (no need to flatten)
        return nested_list[0]
    else:
        # Generate combinations of elements from sublists
        # Apply the specified operation (evaluated as a string) to each combination
        flattened_list = [
            eval(f"{operation.join(map(str, sublist))}")
            for sublist in itertools.product(*nested_list)
        ]
        
        return flattened_list