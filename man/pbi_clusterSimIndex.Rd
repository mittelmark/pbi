\name{pbi_clusterSimIndex} 
\alias{pbi$clusterSimIndex} 
\alias{pbi_clusterSimIndex} 
\title{Compute similarity indices of two clusterings.} 
\description{ 
  This function computes the similarity indices between two clusterings,  
  such as the Rand, Jaccard and Cohen's Kappa index. 
} 
\usage{pbi_clusterSimIndex(v1, v2)} 
\arguments{ 
  \item{v1}{ 
    vector with cluster ids for first clustering 
  } 
  \item{v2}{ 
    vector with cluster ids for second clustering 
  } 
} 
\value{return list with the components ...} 
\examples{ 
  set.seed(123) 
  hcl1=hclust(dist(iris[,1:4])) 
  hcl2=hclust(dist(scale(iris[,1:4]))) 
  round(unlist(pbi_clusterSimIndex(cutree(hcl1,3), 
    cutree(hcl2,3))),2) 
  round(unlist(pbi_clusterSimIndex(cutree(hcl2,3), 
    cutree(hcl1,3))),2) # symmetric 
  round(unlist(pbi_clusterSimIndex(cutree(hcl1,3), 
    sample(cutree(hcl2,3)))),2) # random values I 
  round(unlist(pbi_clusterSimIndex(sample(cutree(hcl1,3)), 
    sample(cutree(hcl2,3)))),2) # random values II 
} 
