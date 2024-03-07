library(purrr)
library(jsonlite)

json_string <- paste(readLines("/workdir/agh227/reference_ECM_signatures/API_Reactome_ECMOrg.txt", warn = FALSE), collapse = "")
# Convert JSON string to a list
json_list <- fromJSON(json_string, simplifyVector = FALSE)

# Extract data using purrr
data <- fromJSON(json_string, flatten = TRUE)

# Print the resulting data frame
print(data)

write.table(data, file = "/workdir/agh227/reference_ECM_signatures/reactome_ECMOrganization.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

