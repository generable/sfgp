% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/models-TSModel.R
\name{TSModel}
\alias{TSModel}
\title{Time series (longitudinal) model class (R6 class)}
\description{
Time series (longitudinal) model class (R6 class)

Time series (longitudinal) model class (R6 class)
}
\section{Super classes}{
\code{\link[sfgp:StanCodeCreator]{sfgp::StanCodeCreator}} -> \code{\link[sfgp:StanModel]{sfgp::StanModel}} -> \code{TSModel}
}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{term_list}}{The additive model terms.}

\item{\code{y_var}}{Name of the y variable.}

\item{\code{id_var}}{Name of the subject identifier variable.}

\item{\code{prior_sigma}}{Prior for the noise parameter.}

\item{\code{sigma_upper}}{Upper bound for the noise parameter.}

\item{\code{sigma_lower}}{Lower for the noise parameter.}

\item{\code{log_y_cap}}{Upper bound on log scale for creating a capped
predicted signal or predicted observations.}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-TSModel-new}{\code{TSModel$new()}}
\item \href{#method-TSModel-stanname_y}{\code{TSModel$stanname_y()}}
\item \href{#method-TSModel-string}{\code{TSModel$string()}}
\item \href{#method-TSModel-get_delta}{\code{TSModel$get_delta()}}
\item \href{#method-TSModel-term_names}{\code{TSModel$term_names()}}
\item \href{#method-TSModel-create_standata}{\code{TSModel$create_standata()}}
\item \href{#method-TSModel-fit}{\code{TSModel$fit()}}
\item \href{#method-TSModel-clone}{\code{TSModel$clone()}}
}
}
\if{html}{\out{
<details><summary>Inherited methods</summary>
<ul>
<li><span class="pkg-link" data-pkg="sfgp" data-topic="StanCodeCreator" data-id="print"><a href='../../sfgp/html/StanCodeCreator.html#method-StanCodeCreator-print'><code>sfgp::StanCodeCreator$print()</code></a></span></li>
<li><span class="pkg-link" data-pkg="sfgp" data-topic="StanCodeCreator" data-id="stancode_data"><a href='../../sfgp/html/StanCodeCreator.html#method-StanCodeCreator-stancode_data'><code>sfgp::StanCodeCreator$stancode_data()</code></a></span></li>
<li><span class="pkg-link" data-pkg="sfgp" data-topic="StanCodeCreator" data-id="stancode_functions"><a href='../../sfgp/html/StanCodeCreator.html#method-StanCodeCreator-stancode_functions'><code>sfgp::StanCodeCreator$stancode_functions()</code></a></span></li>
<li><span class="pkg-link" data-pkg="sfgp" data-topic="StanCodeCreator" data-id="stancode_gq"><a href='../../sfgp/html/StanCodeCreator.html#method-StanCodeCreator-stancode_gq'><code>sfgp::StanCodeCreator$stancode_gq()</code></a></span></li>
<li><span class="pkg-link" data-pkg="sfgp" data-topic="StanCodeCreator" data-id="stancode_model"><a href='../../sfgp/html/StanCodeCreator.html#method-StanCodeCreator-stancode_model'><code>sfgp::StanCodeCreator$stancode_model()</code></a></span></li>
<li><span class="pkg-link" data-pkg="sfgp" data-topic="StanCodeCreator" data-id="stancode_pars"><a href='../../sfgp/html/StanCodeCreator.html#method-StanCodeCreator-stancode_pars'><code>sfgp::StanCodeCreator$stancode_pars()</code></a></span></li>
<li><span class="pkg-link" data-pkg="sfgp" data-topic="StanCodeCreator" data-id="stancode_tdata"><a href='../../sfgp/html/StanCodeCreator.html#method-StanCodeCreator-stancode_tdata'><code>sfgp::StanCodeCreator$stancode_tdata()</code></a></span></li>
<li><span class="pkg-link" data-pkg="sfgp" data-topic="StanCodeCreator" data-id="stancode_tpars"><a href='../../sfgp/html/StanCodeCreator.html#method-StanCodeCreator-stancode_tpars'><code>sfgp::StanCodeCreator$stancode_tpars()</code></a></span></li>
<li><span class="pkg-link" data-pkg="sfgp" data-topic="StanCodeCreator" data-id="stanfiles_functions"><a href='../../sfgp/html/StanCodeCreator.html#method-StanCodeCreator-stanfiles_functions'><code>sfgp::StanCodeCreator$stanfiles_functions()</code></a></span></li>
<li><span class="pkg-link" data-pkg="sfgp" data-topic="StanModel" data-id="compile"><a href='../../sfgp/html/StanModel.html#method-StanModel-compile'><code>sfgp::StanModel$compile()</code></a></span></li>
<li><span class="pkg-link" data-pkg="sfgp" data-topic="StanModel" data-id="create_stancode"><a href='../../sfgp/html/StanModel.html#method-StanModel-create_stancode'><code>sfgp::StanModel$create_stancode()</code></a></span></li>
<li><span class="pkg-link" data-pkg="sfgp" data-topic="StanModel" data-id="create_stanmodel"><a href='../../sfgp/html/StanModel.html#method-StanModel-create_stanmodel'><code>sfgp::StanModel$create_stanmodel()</code></a></span></li>
<li><span class="pkg-link" data-pkg="sfgp" data-topic="StanModel" data-id="get_stanmodel"><a href='../../sfgp/html/StanModel.html#method-StanModel-get_stanmodel'><code>sfgp::StanModel$get_stanmodel()</code></a></span></li>
</ul>
</details>
}}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-TSModel-new"></a>}}
\if{latex}{\out{\hypertarget{method-TSModel-new}{}}}
\subsection{Method \code{new()}}{
Create model
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{TSModel$new(
  formula,
  id_var = "id",
  compile = TRUE,
  delta = 0,
  baseline = NULL,
  prior_baseline = NULL,
  prior_terms = NULL,
  prior_sigma = "normal(0, 2)",
  sigma_upper = 3,
  sigma_lower = 0
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{formula}}{The model formula determining the terms and the y
variable (longitudinal observation).}

\item{\code{id_var}}{Name of the subject identifier variable.}

\item{\code{compile}}{Should the 'Stan' model code be created and compiled.}

\item{\code{delta}}{Offset for log transform (\code{y_log = log(y + delta)}).}

\item{\code{baseline}}{Baseline term definition. Created automatically based
on \code{id_var} if \code{NULL} (default).}

\item{\code{prior_baseline}}{Prior for the baseline term.}

\item{\code{prior_terms}}{A list with names equal to a subset of the
names of the model terms. Can be used to edit priors of term parameters.}

\item{\code{prior_sigma}}{Prior for sigma}

\item{\code{sigma_upper}}{Upper bound for sigma}

\item{\code{sigma_lower}}{Lower bound for sigma}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-TSModel-stanname_y"></a>}}
\if{latex}{\out{\hypertarget{method-TSModel-stanname_y}{}}}
\subsection{Method \code{stanname_y()}}{
Get name of y variable in Stan code.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{TSModel$stanname_y()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-TSModel-string"></a>}}
\if{latex}{\out{\hypertarget{method-TSModel-string}{}}}
\subsection{Method \code{string()}}{
The model description as a string
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{TSModel$string()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-TSModel-get_delta"></a>}}
\if{latex}{\out{\hypertarget{method-TSModel-get_delta}{}}}
\subsection{Method \code{get_delta()}}{
Get value of \code{delta}.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{TSModel$get_delta()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-TSModel-term_names"></a>}}
\if{latex}{\out{\hypertarget{method-TSModel-term_names}{}}}
\subsection{Method \code{term_names()}}{
Get term names in Stan code.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{TSModel$term_names()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
A character vector with length equal to number of terms.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-TSModel-create_standata"></a>}}
\if{latex}{\out{\hypertarget{method-TSModel-create_standata}{}}}
\subsection{Method \code{create_standata()}}{
Create the 'Stan' data list from a data frame. Performs normalization
on continuous variables that are input to GP terms.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{TSModel$create_standata(
  data,
  term_confs = NULL,
  num_bf = NULL,
  scale_bf = NULL,
  skip_transform = NULL,
  prior_only = FALSE,
  set_transforms = TRUE,
  dataname = "LON"
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{data}}{A data frame.}

\item{\code{term_confs}}{A list that specifies configuration of model terms.
If name of any term is not found from the list, \code{$default_conf()}
of that \code{FormulaTerm} is used.}

\item{\code{num_bf}}{If not \code{NULL}, configurations of all
\code{GPTerm}s are updated with this value.}

\item{\code{scale_bf}}{If not \code{NULL}, configurations of all
\code{GPTerm}s are updated with this value.}

\item{\code{skip_transform}}{Term names whose input transform should be
skipped.}

\item{\code{prior_only}}{Sample from prior only?}

\item{\code{set_transforms}}{If data transforms should be set based on the given
\code{data}. This should be \code{TRUE} when fitting a model, and
\code{FALSE} when computing predictions using GQ.}

\item{\code{dataname}}{Name of dataset.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A list.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-TSModel-fit"></a>}}
\if{latex}{\out{\hypertarget{method-TSModel-fit}{}}}
\subsection{Method \code{fit()}}{
Fit the model.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{TSModel$fit(
  data,
  term_confs = NULL,
  num_bf = NULL,
  scale_bf = NULL,
  skip_transform = NULL,
  prior_only = FALSE,
  ...
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{data}}{A data frame.}

\item{\code{term_confs}}{A list that specifies configuration of model terms.
If name of any term is not found from the list, \code{$default_conf()}
of that \code{FormulaTerm} is used.}

\item{\code{num_bf}}{If not \code{NULL}, configurations of all
\code{GPTerm}s are updated with this value.}

\item{\code{scale_bf}}{If not \code{NULL}, configurations of all
\code{GPTerm}s are updated with this value.}

\item{\code{skip_transform}}{Term names whose input transform should be
skipped.}

\item{\code{prior_only}}{Sample from prior only.}

\item{\code{...}}{Arguments passed to \code{sample} method of the
'CmdStanR' model.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
An \code{\link{TSModelFit}} object.
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-TSModel-clone"></a>}}
\if{latex}{\out{\hypertarget{method-TSModel-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{TSModel$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
