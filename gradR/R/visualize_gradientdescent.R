#' Launch the Optimization Diagnostic Dashboard
#'
#' @description
#' This function launches an interactive Shiny application designed to visualize
#' mathematical optimization. It evaluates first-order gradient descent,
#' second-order Hessian conditions (saddle point diagnostics), and simulates
#' Stochastic Gradient Descent (SGD) with learning rate decay.
#'
#' @return Opens a local Shiny application in your default web browser.
#' @export
#'
#' @import shiny
#' @import ggplot2
#' @import plotly
#' @import dplyr
#'
#' @examples
#' \dontrun{
#' # Launch the interactive diagnostic tool
#' visualize_gradientdescent()
#' }
visualize_gradientdescent <- function() {

  # Ensure libraries are loaded when the function runs from the package
  library(shiny)
  library(ggplot2)
  library(plotly)
  library(dplyr)

  # ---------------------------------------------------------
  # UI DEFINITION
  # ---------------------------------------------------------
  ui <- fluidPage(
    withMathJax(),

    titlePanel("Theoretical Gradient Descent & Auto-Tuner"),

    sidebarLayout(
      sidebarPanel(
        h4("1. Define the Landscape"),
        textInput("equation", "Objective Function f(x,y):", value = "x^2 + 2*y^2"),
        helpText("Try the Saddle Point: 'x^2 - y^2', or Convex: 'x^2 + 2*y^2'"),

        hr(),
        h4("2. Optimization Mode"),
        radioButtons("tuning_mode", "Select Mode:",
                     choices = c("Manual Tuning", "Grid Search (Auto-Tune)"),
                     selected = "Manual Tuning"),

        conditionalPanel(
          condition = "input.tuning_mode == 'Manual Tuning'",
          sliderInput("learning_rate", "Initial Learning Rate (\\(\\alpha_0\\)):", min = 0.001, max = 0.5, value = 0.05, step = 0.005)
        ),
        conditionalPanel(
          condition = "input.tuning_mode == 'Grid Search (Auto-Tune)'",
          helpText("Simulating 5 parallel paths with \\(\\alpha_0\\) values: 0.01, 0.05, 0.1, 0.2, 0.4")
        ),

        hr(),
        h4("3. Stochastic Sampling (Noise)"),
        checkboxInput("add_noise", "Enable Stochastic Sampling", value = FALSE),
        conditionalPanel(
          condition = "input.add_noise == true",
          sliderInput("noise_sd", "Sampling Standard Deviation (\\(\\sigma\\)):", min = 0.01, max = 1.0, value = 0.2, step = 0.01),
          sliderInput("decay_rate", "Decay Rate (\\(k\\)):", min = 0, max = 0.2, value = 0.05, step = 0.01),
          helpText("Simulates SGD. Decay shrinks the learning rate over time so the algorithm can settle despite the noise.")
        ),

        hr(),
        h4("4. Algorithm Parameters"),
        sliderInput("start_x", "Starting X Position:", min = -5, max = 5, value = -4, step = 0.5),
        sliderInput("start_y", "Starting Y Position:", min = -5, max = 5, value = 4, step = 0.5),
        sliderInput("iterations", "Max Iterations:", min = 1, max = 200, value = 50, step = 1),
        numericInput("tolerance", "Convergence Tolerance (\\(\\epsilon\\)):", value = 0.01, step = 0.005),

        hr(),
        h4("5. Visualization"),
        radioButtons("view_mode", "Landscape View:", choices = c("2D Contour Map", "3D Surface Plot"), selected = "2D Contour Map"),
        conditionalPanel(
          condition = "input.view_mode == '2D Contour Map'",
          checkboxInput("show_vectors", "Overlay Gradient Vector Field", value = TRUE)
        )
      ),

      mainPanel(
        tags$style(type="text/css", ".shiny-output-error { color: red; font-weight: bold;}"),

        wellPanel(
          h4("Algorithm Diagnostics"),
          uiOutput("math_readout"),
          hr(),
          h4("Second-Order Optimality Test"),
          uiOutput("hessian_readout")
        ),

        h4("Loss Landscape & Descent Paths"),
        plotlyOutput("landscapePlot", height = "450px"),

        hr(),
        h4("Comparative Loss Curve"),
        plotOutput("lossCurve", height = "250px")
      )
    )
  )

  # ---------------------------------------------------------
  # SERVER LOGIC
  # ---------------------------------------------------------
  server <- function(input, output) {

    parsed_function <- reactive({
      req(input$equation)
      tryCatch({
        deriv(parse(text = input$equation), c("x", "y"), function.arg = TRUE, hessian = TRUE)
      }, error = function(e) { return(NULL) })
    })

    gd_engine <- reactive({
      f_grad <- parsed_function()
      validate(need(!is.null(f_grad), "Error: Invalid mathematical expression."))

      start_x <- input$start_x
      start_y <- input$start_y
      epsilon <- input$tolerance
      max_iter <- input$iterations

      noise_sd <- ifelse(input$add_noise, input$noise_sd, 0)
      decay_rate <- ifelse(input$add_noise, input$decay_rate, 0)

      alphas <- if (input$tuning_mode == "Manual Tuning") {
        input$learning_rate
      } else {
        c(0.01, 0.05, 0.1, 0.2, 0.4)
      }

      all_histories <- data.frame()
      summary_list <- list()
      initial_grad <- c(0, 0)

      for (alpha in alphas) {
        current_x <- start_x
        current_y <- start_y

        history <- data.frame(
          iteration = 0, x = current_x, y = current_y,
          z = suppressWarnings(as.numeric(f_grad(current_x, current_y))),
          grad_mag = NA, alpha = as.factor(alpha)
        )

        converged <- FALSE
        status_msg <- paste("Reached Max Iterations (", max_iter, ")", sep="")
        iter_taken <- max_iter

        for (i in 1:max_iter) {
          eval_result <- suppressWarnings(f_grad(current_x, current_y))
          grad_vector <- attr(eval_result, "gradient")[1, ]

          if (alpha == alphas[1] && i == 1) initial_grad <- grad_vector

          if (input$add_noise) {
            grad_vector["x"] <- grad_vector["x"] + rnorm(1, mean = 0, sd = noise_sd)
            grad_vector["y"] <- grad_vector["y"] + rnorm(1, mean = 0, sd = noise_sd)
          }

          if (is.na(grad_vector["x"]) || is.infinite(grad_vector["x"])) {
            status_msg <- "Diverged (Gradient Exploded!)"
            iter_taken <- i
            break
          }

          grad_magnitude <- sqrt(grad_vector["x"]^2 + grad_vector["y"]^2)
          history$grad_mag[nrow(history)] <- grad_magnitude

          if (!is.na(grad_magnitude) && grad_magnitude < epsilon && !input$add_noise) {
            status_msg <- paste("Converged in", i, "steps")
            converged <- TRUE
            iter_taken <- i
            break
          }

          current_alpha <- alpha / (1 + decay_rate * i)

          current_x <- current_x - (current_alpha * grad_vector["x"])
          current_y <- current_y - (current_alpha * grad_vector["y"])

          new_z <- suppressWarnings(as.numeric(f_grad(current_x, current_y)))
          history <- rbind(history, data.frame(iteration = i, x = current_x, y = current_y, z = new_z, grad_mag = NA, alpha = as.factor(alpha)))
        }

        final_eval <- suppressWarnings(f_grad(current_x, current_y))
        final_hessian <- attr(final_eval, "hessian")[1, , ]
        eigenvals <- eigen(final_hessian)$values

        if (all(eigenvals > 1e-5)) {
          point_type <- "Strict Local Minimum (Positive Definite)"
          color_code <- "green"
        } else if (all(eigenvals < -1e-5)) {
          point_type <- "Strict Local Maximum (Negative Definite)"
          color_code <- "red"
        } else if (eigenvals[1] * eigenvals[2] < 0) {
          point_type <- "Saddle Point (Indefinite Matrix)"
          color_code <- "orange"
        } else {
          point_type <- "Inconclusive (Semi-Definite)"
          color_code <- "gray"
        }

        all_histories <- rbind(all_histories, history)
        summary_list[[as.character(alpha)]] <- list(
          alpha = alpha, converged = converged, iters = iter_taken,
          status = status_msg, hessian = final_hessian,
          eigenvals = eigenvals, point_type = point_type, color = color_code,
          final_x = current_x, final_y = current_y, final_z = new_z
        )
      }

      best_summary <- summary_list[[1]]
      if (length(alphas) > 1) {
        converged_runs <- Filter(function(x) x$converged, summary_list)
        if (length(converged_runs) > 0) {
          best_summary <- converged_runs[[which.min(sapply(converged_runs, function(x) x$iters))]]
        } else {
          best_summary$status <- "None converged. Try different parameters."
        }
      }

      return(list(history = all_histories, best_summary = best_summary,
                  init_grad = initial_grad, is_grid = (input$tuning_mode == "Grid Search (Auto-Tune)"),
                  add_noise = input$add_noise))
    })

    output$math_readout <- renderUI({
      engine <- gd_engine()
      res <- engine$best_summary
      init_g <- engine$init_grad

      mode_text <- if (engine$is_grid) {
        paste0("<p><b style='color:#007bc2; font-size:16px;'>Grid Search Result:</b> Optimal \\(\\alpha_0\\) = <b>", res$alpha, "</b> (", res$status, ")</p>")
      } else {
        paste0("<p><b>Status:</b> ", res$status, "</p>")
      }

      formula_text <- if (engine$add_noise) {
        "$$\\mathbf{x}_{i+1} = \\mathbf{x}_i - \\left( \\frac{\\alpha_0}{1 + k \\cdot i} \\right) \\left( \\nabla f(\\mathbf{x}_i) + \\mathcal{N}(0, \\sigma^2) \\right)$$"
      } else {
        "$$\\mathbf{x}_{i+1} = \\mathbf{x}_i - \\alpha \\nabla f(\\mathbf{x}_i)$$"
      }

      withMathJax(HTML(paste0(
        mode_text,
        "<p><b>Theoretical Update Rule:</b> ", formula_text, "</p>",
        "<p><b>Step 1 Calculus:</b> Initial Gradient Vector \\(\\nabla f(x_0, y_0)\\) = [", round(init_g["x"], 3), ", ", round(init_g["y"], 3), "]</p>",
        "<p><b>Stopping Criterion:</b> Break when \\(||\\nabla f|| < \\epsilon\\)</p>",
        "<hr>",
        "<p><b>Final Resting Point (x, y):</b> [", round(res$final_x, 4), ", ", round(res$final_y, 4), "]</p>",
        "<p><b style='color:#d9534f; font-size:16px;'>Minimum Function Value f(x,y):</b> <b style='font-size:16px; color:#d9534f;'>", round(res$final_z, 6), "</b></p>"
      )))
    })

    output$hessian_readout <- renderUI({
      res <- gd_engine()$best_summary
      H <- res$hessian
      eig <- res$eigenvals

      h_matrix_tex <- paste0(
        "$$ H = \\begin{bmatrix} ", round(H[1,1], 3), " & ", round(H[1,2], 3), " \\\\ ",
        round(H[2,1], 3), " & ", round(H[2,2], 3), " \\end{bmatrix} $$"
      )

      withMathJax(HTML(paste0(
        "<p><b>Final Hessian Matrix:</b> Evaluates second-order curvature at the resting point:</p>",
        h_matrix_tex,
        "<p><b>Eigenvalues (\\(\\lambda_1, \\lambda_2\\)):</b> [", round(eig[1], 3), ", ", round(eig[2], 3), "]</p>",
        "<p style='font-size:16px;'><b>Topological Classification:</b> <b style='color:", res$color, "'>", res$point_type, "</b></p>"
      )))
    })

    output$landscapePlot <- renderPlotly({
      f_grad <- parsed_function()
      history <- gd_engine()$history

      grid_limit <- min(max(abs(c(history$x, history$y)), na.rm = TRUE) + 2, 15)
      x_seq <- seq(-grid_limit, grid_limit, length.out = 60)
      y_seq <- seq(-grid_limit, grid_limit, length.out = 60)

      z_matrix <- matrix(nrow = length(x_seq), ncol = length(y_seq))
      for (i in 1:length(x_seq)) {
        for (j in 1:length(y_seq)) {
          z_matrix[i, j] <- suppressWarnings(as.numeric(f_grad(x_seq[i], y_seq[j])))
        }
      }

      if (input$view_mode == "2D Contour Map") {
        p <- plot_ly(x = ~x_seq, y = ~y_seq, z = t(z_matrix), type = "contour",
                     colorscale = "Viridis", contours = list(showlabels = TRUE))

        if (input$show_vectors) {
          vec_seq <- seq(-grid_limit, grid_limit, length.out = 15)
          vec_grid <- expand.grid(x = vec_seq, y = vec_seq)
          arrows <- lapply(1:nrow(vec_grid), function(i) {
            eval_res <- suppressWarnings(f_grad(vec_grid$x[i], vec_grid$y[i]))
            g_vec <- attr(eval_res, "gradient")[1, ]
            mag <- sqrt(g_vec["x"]^2 + g_vec["y"]^2)
            if(is.na(mag) || mag == 0) return(NULL)
            scale_factor <- (grid_limit * 2 / 15) * 0.4 / mag
            list(x = vec_grid$x[i] + (g_vec["x"] * scale_factor), y = vec_grid$y[i] + (g_vec["y"] * scale_factor),
                 ax = vec_grid$x[i], ay = vec_grid$y[i], xref = "x", yref = "y", axref = "x", ayref = "y",
                 showarrow = TRUE, arrowhead = 2, arrowsize = 0.8, arrowcolor = "rgba(255,255,255,0.6)")
          })
          p <- p %>% layout(annotations = Filter(Negate(is.null), arrows))
        }

        p %>% add_trace(
          data = history, x = ~x, y = ~y, color = ~alpha, type = "scatter", mode = "lines+markers",
          line = list(width = 3), marker = list(size = 5), text = ~paste("Iteration:", iteration)
        ) %>% layout(xaxis = list(title = "X"), yaxis = list(title = "Y"))

      } else {
        plot_ly(x = ~x_seq, y = ~y_seq, z = t(z_matrix)) %>%
          add_surface(colorscale = "Viridis", opacity = 0.8) %>%
          add_trace(
            data = history, x = ~x, y = ~y, z = ~z, color = ~alpha, type = "scatter3d", mode = "lines+markers",
            line = list(width = 5), marker = list(size = 4)
          )
      }
    })

    output$lossCurve <- renderPlot({
      history <- gd_engine()$history
      safe_history <- history %>% filter(is.finite(z))

      ggplot(safe_history, aes(x = iteration, y = z, color = alpha, group = alpha)) +
        geom_line(linewidth = 1) +
        geom_point(size = 2.5) +
        theme_minimal(base_size = 14) +
        labs(x = "Iteration Step", y = "Loss Value: f(x,y)", color = "Learning Rate (α0)") +
        theme(legend.position = "right", panel.grid.minor = element_blank())
    })
  }

  # Run the application
  shinyApp(ui = ui, server = server)
}
