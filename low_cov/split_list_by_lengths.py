def split_list_by_lengths(input_list, lengths_list):
    """
    Split a list into sublists of specified lengths.
    
    Args:
        input_list (list): The list to be split.
        lengths_list (list): A list of integers representing the lengths of the sublists.
    
    Returns:
        list: A list of sublists created based on the specified lengths.
    """
    split_list = []  # Initialize an empty list to store the sublists
    start = 0  # Initialize the starting index
    
    # Iterate through the lengths in the lengths_list
    for length in lengths_list:
        # Extract a sublist from input_list
        sublist = input_list[start:start + length]
        
        # Append the sublist to the result list
        split_list.append(sublist)
        
        # Update the starting index for the next iteration
        start += length
    
    return split_list