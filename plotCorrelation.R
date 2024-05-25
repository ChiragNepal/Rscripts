# Set seed for reproducibility
set.seed(123)

# Generate a data matrix with 100 observations of 5 variables
data_matrix <- matrix(rnorm(500), nrow=100, ncol=5)

# Convert the matrix to a data frame for easier handling
data_df <- as.data.frame(data_matrix)

# Give meaningful names to the columns
colnames(data_df) <- c("Var1", "Var2", "Var3", "Var4", "Var5")

# Display the first few rows of the data frame
head(data_df)

# Calculate the correlation matrix
cor_matrix <- cor(data_df)

# Display the correlation matrix
print(cor_matrix)

# Install and load the corrplot package
install.packages("corrplot")
library(corrplot)

# Create the correlation plot
corrplot(cor_matrix, method = "circle", type = "upper", tl.col = "black", tl.srt = 45, 
         addCoef.col = "black", number.cex = 0.7)


