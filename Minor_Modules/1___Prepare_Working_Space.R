### ---------------------
### Prepare the workspace
### ---------------------

# Clear the R workspace
rm(list = ls(all=T))

# Allocate larger than usual memory and 
memory.size(5000000)
options(java.parameters = "-Xmx12288m")