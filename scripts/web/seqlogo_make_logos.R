require(ggplot2)
require(ggseqlogo)

# Get the input pwm and the output file name
args = commandArgs(trailingOnly=TRUE)
input_pwm = args[1]
output_file = args[2]

# Read a .csv file containing the nucleotide frequencies
data = read.csv(input_pwm, header=FALSE)
# Transverse the nucleotide frequencies table
data = t(data)
# Assign nucleotides to each row
rownames(data) <- c("A", "C", "G", "T")
# Make the plot 
ylim = c(0,2)
pwm_logo = ggseqlogo(data)
# Save the plot
len_data = (length(data)*0.15)
ggsave(filename=output_file, plot=pwm_logo, height=2, width=len_data, limitsize=FALSE)
#ggsave(filename=output_file, plot=pwm_logo, limitsize=TRUE)