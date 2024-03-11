\name{pbi_clusterSilhouette} 
\alias{pbi$clusterSilhouette} 
\alias{pbi_clusterSilhouette} 
\title{Determine Cluster Silhouette Index} 
\description{ 
    This function can be used to determine the strength of a given clustering. 
    A silhouette index of larger than 0.7 indicates a strong,  
    values between 0.5 and 0.7 a reasonable, between 0.25 and 0.5 a weak and below 
    0.25 no reasonable structure. 
} 
\usage{pbi_clusterSilhouette(x,D)} 
\arguments{ 
  \item{x}{cluster ids} 
  \item{D}{a distance matrix object} 
} 
\value{return a list with the components avg.width (the silhouette index value for the complete) and clus.avg.widths the indices for the individual clusters.} 
\examples{ 
  D=dist(scale(iris[,1:4])) 
  hcl=hclust(D) 
  pbi_clusterSilhouette(cutree(hcl,2),D) 
  pbi_clusterSilhouette(cutree(hcl,3),D)$avg.width 
  pbi_clusterSilhouette(cutree(hcl,4),D)$avg.width 
  pbi_clusterSilhouette(cutree(hcl,5),D)$avg.width 
  # three cluster seems to be the best 
} 
