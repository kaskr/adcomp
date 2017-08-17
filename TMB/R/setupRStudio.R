## Experimental RStudio integration
setupRStudio <- function(file="~/.Rprofile",
                         snipRfile="~/.R/snippets/r.snippets",
                         snipCppfile="~/.R/snippets/c_cpp.snippets") {
code <- 
'
######## TMB - setup RStudio
setHook(packageEvent("TMB", "onLoad"),
        function(...) {
            if("tools:rstudio" %in% search()) {
                tmb.env <- asNamespace("TMB")
                compile.orig <- tmb.env$compile
                unlockBinding("compile", tmb.env)
                ## Rstudio handle compilation errors:
                rs.env <- as.environment("tools:rstudio")
                tmb.env$compile <- function(file,...) {
                    .Call("rs_sourceCppOnBuild",
                          file, FALSE, FALSE)
                    status <- try( compile.orig(file, ...) )
                    succeeded <- (status == 0)
                    .Call("rs_sourceCppOnBuildComplete",
                          succeeded, "")
                    if(!succeeded) stop("Compilation failed")
                    status
                }
                ## Bind "sourceCpp" button to TMB compile
                rcpp.env <- asNamespace("Rcpp")
                unlockBinding("sourceCpp", rcpp.env)
                rcpp.env$sourceCpp <- tmb.env$compile
            }
        } )
'
mess <- c("You are about to setup Rstudio with TMB.",
          "Changes will be added to the file:",
          "",
          file,
          "")
invisible(lapply(mess, cat, "\n"))
ans <- readline("OK? (yes/no) ")
if(ans == "yes") {
    chk <- grepl(
        "######## TMB - setup RStudio",
        readLines(file)
    )
    if(any(chk))
        message("Skipping because changes seem to have been made already.")
    else {
        cat(code, file=file, append=TRUE)
        message("Please re-start RStudio for the changes to take place.")
    }
} else {
    message("Dropping out")
}

## Experimental RStudio TMB snippet integration
## Gavin Fay & Andrea Havron
## 
  rsnips <- getRsnips()
  headers <- grep("snippet",rsnips)
  rsnips[-headers] <- paste("\t ",rsnips[-headers],sep="")
  cppsnips <- getCppsnips()
  cheaders <- grep("snippet",cppsnips)
  cppsnips[-cheaders] <- paste("\t ",cppsnips[-cheaders],sep="")
  
  mess <- c("You are about to setup snippets for TMB.",
            "Changes will be added to the files:",
            "",
            snipRfile,
            "",
            snipCppfile,
            "")
  invisible(lapply(mess, cat, "\n"))
  ans <- readline("OK? (yes/no) ")
  if(ans == "yes") {
    chk <- grepl(
      rsnips[headers][1],
      readLines(snipRfile)
    )
    if(any(chk))
      message("Skipping because changes seem to have been made already.")
    else {
      #cat(snips, file=outfile, append=TRUE)
      write.table(rsnips,file=snipRfile,quote=FALSE,
                  row.names=FALSE,col.names=FALSE,append=TRUE)
      write.table(cppsnips,file=snipCppfile,quote=FALSE,
                  row.names=FALSE,col.names=FALSE,append=TRUE)
      message("Please re-start RStudio for the changes to take place.")
    }
  } else {
    message("Dropping out")
  }
  
  invisible(NULL)
}



#R snippets, paste as text, replace \$ with \\$
getRsnips <- function() {
  snips <- 
'
snippet tmb.sp 
library(TMB)
library(RandomFields)
# Compile the C++ file
compile(${1:"spatial_spde.cpp"})
# Dynamically link the C++ code
dyn.load(dynlib(${2:"spatial_spde"}))     
  
loc <- ${3:matrix(runif(200), nrow = 100, ncol = 2)}
  
#Simulate Field
RMmodel <- RMgauss(var=1, scale=0.5) 
omega <- array(RFsimulate(model=RMmodel, x=loc[,1], y=loc[,2])@data[,1], dim=c(100,100))
  
#Create mesh
mesh <- inla.mesh.2d(loc, max.edge = ${4:c(0.2,0.5)}, ${5:...})
spde <- inla.spde2.matern(mesh, alpha=2)
  
#Create Data and Parameter Lists
data <- list(${6: y = rpois(100, exp(10 + omega))}, v_i = mesh\\$idx\\$loc-1)
data\\$spde <- spde\\$param.inla[c("M0","M1","M2")]
params <- list(${7: eta = 0}, ln_kappa = 0, ln_tau = 0, ${8:Omega} = rep(0, spde\\$n.spde))
random <- ${9:"Omega"}
  
#make the function
f <- MakeADFun(data = data, parameters = params, random = random)
#optimize parameters
fit <- nlminb( f\\$par , f\\$fn , f\\$gr  )
#generate the report
rep <- f\\$rep()
  
#get standard errors
sdrep <- sdreport(f)
  
#check convergence
sdrep\\$pdHess #TRUE
f\\$gr(fit\\$par) #values close to 0 or less
  
snippet tmb.mf
MakeADFun(${1:list}, ${2:list}, random = ${3:vector}, ${4:...})

snippet tmb.nlm
nlminb(${1:par}, ${2:fn}, ${3:gr}, ${4:...})

snippet tmb.load
dyn.load(dynlib(${1:"dll"}))
  
snippet tmb.unload
dyn.unload(dynlib(${1:"dll"}))
  
snippet tmb.spde
inla.spde2.matern(${1:mesh}, alpha = ${2:2}, ${3:...})
  
snippet tmb.mesh
inla.mesh.2d(loc = ${1:matrix}, offset = ${2:value or vector}, max.edge = ${3:value or vector}, cutoff = ${4:value}, ${5:...})
  

snippet tmb.rscript
library(TMB)    

# Compile the C++ file
compile(${1:"tutorial.cpp"})

# Dynamically link the C++ code
dyn.load(dynlib(${2:"tutorial"}))               

#set up data object
data <- list(${3:x = rnorm(100,0,10)})

#define parameters of the model
params <- list(${4:mu=0}, ${5:sigma=1})

#specify the parameters that are random effects
#if no random effects, specify NULL
random <- ${6:NULL}

#make the function
f <- MakeADFun(data = data, parameters = params, random = random)

#optimize parameters
fit <- nlminb( f\\$par , f\\$fn , f\\$gr , lower = ${7:c(-10,0.0)}, upper = ${8:c(10.0,10.0)} )

#generate the report
rep <- f\\$rep()
#get standard errors
sdrep <- sdreport(f)

#check convergence
sdrep\\$pdHess #TRUE
f\\$gr(fit\\$par) #values close to 0
'
  # Following lines need improving to avoid writing to temporary file.
  write(snips,file="snips.out")  
  snips <- readLines("snips.out")
  file.remove("snips.out")
  return(snips)
}


#R snippets, paste as text, replace \$ with \\$
getCppsnips <- function() {
  snips <- 
'
snippet spatial_spde.cpp
  
#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
using namespace R_inla;
using namespace density;
using namespace Eigen;

// Data
DATA_VECTOR( y );        // Number of factors we want to estimate -
DATA_IVECTOR( v_i ); // Convert site s to vertex x using SPDE method
DATA_STRUCT(spde,spde_t);


PARAMETER( eta );
PARAMETER( ln_kappa );
PARAMETER( ln_tau );

// Random effects
PARAMETER_VECTOR( Omega );

int n_y = y.size();
vector<Type> nll_comp(2);
nll_comp.setZero();

// Derived quantities
Type range = sqrt(8) / exp( ln_kappa ); //distance at which correlation is ~10%
Type sigma2 = exp( - log( 4 * M_PI) - (2*ln_kappa) - (2*ln_tau) ); //marginal spatial variance

SparseMatrix<Type> Q = Q_spde(spde,exp(ln_kappa));
nll_comp(0) += SCALE( GMRF(Q), 1/exp(ln_tau) )( Omega );


vector<Type> ln_mean(n_y);
for(int i=0; i<n_y; i++){
ln_mean(i) = eta + Omega(v_i(i));
nll_comp(1) -= dpois(y(i), exp(ln_mean(i)), true);
}
Type nll = nll_comp.sum();

REPORT( nll_comp );
REPORT( nll );
REPORT( range );
REPORT( sigma2 );

return nll;
}

'
  # Following lines need improving to avoid writing to temporary file.
  write(snips,file="snips.out")  
  snips <- readLines("snips.out")
  file.remove("snips.out")
  return(snips)
}
