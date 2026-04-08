def _validate_site_file_names(site_file_name, file_names):
    """
    This checks that the file names are all unique or all the same.
    """
    unique_names = list(dict.fromkeys(file_names))
    if 1 < len(unique_names) < len(file_names):
        msg = (
            "There are some duplicate file names in {0}. Use either the same file name "
            "for all sites or a different file name for each site."
        ).format(site_file_name)
        raise ValueError(msg)
