extract_binding_energy <- function(file) {
    require(stringr)

    lines <- readLines(file)
    # Find the first REMARK VINA RESULT line
    result_line <- grep("REMARK VINA RESULT", lines, value = TRUE)[1]
    if (!is.null(result_line)) {
        # Extract the binding energy (first number after VINA RESULT)
        energy <- as.numeric(str_extract(result_line, "-?\\d+\\.\\d+"))
        return(energy)
    }
    return(NA)
}