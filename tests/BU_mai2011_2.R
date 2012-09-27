library("planor")
robpl1r4.fac <- planor.factors(factors=list(
                                 row1=1:2,
                                 row2=1:2,
                                 col1=1:2,
                                 col2=1:2,
                                 nsoil=c("curd","Saint-Paulin"),
                                 qsoil=c("10 mg","100 mg"),
                                 cbact=c("3%","6%"),
                                 Tact=c("15 mn","30 mn"),
                                 conc=c("1%","3%"),
                                 brush=c("strong","weak"),
                                 rough=c(0.25,0.75),
                                 nat=1:2),
                               block=~row1+row2+col1+col2,
                               hierarchy=list(~nsoil/(col1*col2*row2),
                                 ~cbact/(col1*col2*row2),
                                 ~Tact/(row1*row2),
                                 ~conc/(row1*row2),
                                 ~brush/(row1*row2*col1)))

robpl1r4.mod <- planor.model( listofmodels=list(
  c( ~row2 + (nsoil+qsoil+cbact+Tact+conc+brush+rough+nat)^2, ~nsoil+qsoil+cbact+Tact+conc+brush+rough+nat), 
  c( ~col1*col2*row2 + row1*row2*col1,                        ~rough+nat)) )

robpl1r4.key <- planor.designkey(factors=robpl1r4.fac, 
                                 model=robpl1r4.mod, 
                                 nunits=16, 
                                 base=~row1+row2+col1+col2 )

summary(robpl1r4.key)
