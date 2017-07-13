\name{madras}
\alias{madras}
\docType{data}
\title{
Madras Longitudinal Schizophrenia Study: Thought Disorder Subset
}
\description{
\code{madras} contains a subset of the data from the Madras Longitudinal Schizophrenia Study, which collected monthly symptom data on 86 schizophrenia patients after their initial hospitalization. The primary question of interest is whether subjects with an older age-at-onset tend to recover more or less quickly, and whether female patients recover more or less quickly. Recovery is measured by a reduction in the presentation of symptoms.
}
\format{
  A data frame with 922 observations on the following 5 variables.
  \describe{
    \item{\code{thought}}{an indicator for thought disorders}
    \item{\code{age}}{an indicator for age-at-onset \eqn{>= 20} years}
    \item{\code{gender}}{an indicator for female gender}
    \item{\code{month}}{months since hospitalization}
    \item{\code{id}}{a unique patient identifier}
  }
}

\examples{
\dontrun{
# Replication of Table 4 from Schildcrout and Heagerty 2007.
data(madras)
madras2 = madras[,c('id','thought','month')]
madras2$gender = (madras$gender==0)*1 # Male
madras2$age = (madras$age==0)*1       # 1= <20

mLV        = mm(thought~month + age + gender + age:month + gender:month,
                lv.formula=~1, data=madras2,id=id)
mT         = mm(thought~month + age + gender + age:month + gender:month, 
                t.formula=~1, data=madras2,id=id)
mTLV       = mm(thought~month + age + gender + age:month + gender:month, 
                lv.formula=~1, t.formula=~1, data=madras2,id=id)
mTLV.lv.g  = mm(thought~month + age + gender + age:month + gender:month, 
                lv.formula=~gender, t.formula=~1, data=madras2,id=id)
mTLV.t.g   = mm(thought~month + age + gender + age:month + gender:month, 
                lv.formula=~1, t.formula=~gender, data=madras2,id=id)
mTLV.lvt.g = mm(thought~month + age + gender + age:month + gender:month,
                t.formula=~gender,lv.formula=~gender, 
                data=madras2,id=id)
}
}

\source{
ADD SOURCE
}


