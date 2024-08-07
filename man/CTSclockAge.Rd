% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CTSclockAge.R
\name{CTSclockAge}
\alias{CTSclockAge}
\title{Predict DNA methylation age using CTS clocks}
\usage{
CTSclockAge(
  data.m,
  CTSclock = c("Neu-In", "Neu-Sin", "Glia-In", "Glia-Sin", "Brain", "Hep", "Liver"),
  dataType = c("bulk", "sorted"),
  CTF.m = NULL,
  tissue = c("brain", "otherTissue"),
  coreNum = NULL
)
}
\arguments{
\item{data.m}{A DNAm matrix (row: CpGs, column: samples) of the samples you want to get
DNAm age predicted by a CTS clock.}

\item{CTSclock}{Which CTSclock you want to use to predict DNAm age ('Neu-In', 'Neu-Sin',
'Glia-In', 'Glia-Sin', 'Brain', 'Hep', 'Liver').}

\item{dataType}{Type of the samples in your DNAm data ('bulk' or 'sorted').}

\item{CTF.m}{Cell type fraction matrix (row: samples in the same order with data.m,
column: cell types).}

\item{tissue}{What tissue are your samples from ('brain', 'otherTissue').}

\item{coreNum}{The number of cores to use. If set to NULL then it will be the number
generated by ceiling(parallel::detectCores()/2).}
}
\value{
A vecter of the predicted DNAm ages.
}
\description{
This is the function for computing DNA methylation age using CTS (cell type
specific) clocks. The inputs including DNAm matrix, the CTS clock you want to
use, is your DNAm data from bulk tissue sample or sorted cell sample, cell
type fraction matrix if you want to use Neu-In/Glia-In/Brain clock, tissue of
your DNAm data samples and the number of cores if you want to do parallel computing.
}
\examples{
data(MurphyGSE88890)
agePred.v = CTSclockAge(beta.m, CTSclock = 'Neu-In', dataType = 'bulk',
CTF.m = NULL, tissue = 'brain')
plot(phenotype.df$Age, agePred.v)

data(ExampleData_Liver)
agePred.v = CTSclockAge(Test.m, CTSclock = 'Hep', dataType = 'bulk',
CTF.m = NULL, tissue = 'otherTissue')
plot(Age, agePred.v)


}
\references{
Huige Tong, Xiaolong Guo, Qi Luo and Andrew E Teschendorff. 2024. "Cell-type specific epigenetic clocks to quantify biological age at cell-type resolution." Submitted.
}
