library(BIGL)
library(rgl)
library(DT)
set.seed(314159265)

shinyServer(function(input, output, session) {

  doseGrid <- reactive({

    if (input$logScale) {
      dose1 <- dose2 <- round(c(0, 3^(-6:0)), 4)
    } else {
      dose1 <- dose2 <- round(seq(0, 3, length.out = 7), 4)
    }

    doseGrid <- expand.grid(list("d1" = sort(dose1),
                                 "d2" = sort(dose2)))
  })

  pars <- reactive({
    c("h1" = input$h1,
      "h2" = input$h2,
      "b" = input$b,
      "m1" = input$m1,
      "m2" = input$m2,
      "e1" = log(input$e1),
      "e2" = log(input$e2))
  })

  ## Default parameters
  observeEvent(input$defaultpars, {
    updateNumericInput(session, "h1", value = 1)
    updateNumericInput(session, "h2", value = 1)
    updateNumericInput(session, "b", value = 0)
    updateNumericInput(session, "m1", value = 1)
    updateNumericInput(session, "m2", value = 1)
    updateNumericInput(session, "e1", value = 0.1)
    updateNumericInput(session, "e2", value = 0.1)
    updateNumericInput(session, "noise", value = 0)
  })

  ## Agonist and partial agonist
  observeEvent(input$pagonist, {
    updateNumericInput(session, "h1", value = 1)
    updateNumericInput(session, "h2", value = 1)
    updateNumericInput(session, "b", value = 0)
    updateNumericInput(session, "m1", value = 1)
    updateNumericInput(session, "m2", value = 0.5)
    updateNumericInput(session, "e1", value = 0.1)
    updateNumericInput(session, "e2", value = 0.1)
    updateNumericInput(session, "noise", value = 0)
  })

  ## Agonist and antagonist
  observeEvent(input$antagonist, {
    updateNumericInput(session, "h1", value = 1)
    updateNumericInput(session, "h2", value = 1)
    updateNumericInput(session, "b", value = 1)
    updateNumericInput(session, "m1", value = 2)
    updateNumericInput(session, "m2", value = 0)
    updateNumericInput(session, "e1", value = 0.1)
    updateNumericInput(session, "e2", value = 0.1)
    updateNumericInput(session, "noise", value = 0.0001)
  })

  comp <- reactive({
    a <- BIGL:::generalizedLoewe(doseGrid(), pars())
    trueo <- with(a$occupancy, (d1/exp(pars()["e1"]) + d2/exp(pars()["e2"]))^pars()["h1"] /
                               (1 + (d1/exp(pars()["e1"]) + d2/exp(pars()["e2"]))^pars()["h1"]))
    a$occupancy$occupancy <- trueo
    a
  })

  dr <- reactive({
    data <- cbind(comp()$occupancy[, c("d1", "d2")], "effect" = comp()$response[-1])
    data <- rbind(c(0, 0, pars()["b"]), data)

    data <- data[rep(row.names(data), input$replicates), ]
    data$effect <- data$effect + input$noise * rnorm(nrow(data), 0, 1)
    data
  })

  fit <- reactive({

    if (input$null == "stdloewe") {
      constraints <- list("matrix" = c(0, 0, 0, 1, -1, 0, 0),
                          "vector" = 0)
    } else {
      constraints <- NULL
    }

    fitMarginals(dr(), method = "nlslm",
                 control = list(maxiter = 200),
                 constraints = constraints)
  })

  compE <- reactive({
    BIGL:::generalizedLoewe(doseGrid(), fit()$coef)
  })

  ## Coefficient table
  output$coefs <- renderTable({
    coefs <- t(rbind(pars()[c("h1", "h2", "e1", "e2", "b", "m1", "m2")],
                     fit()$coef[c("h1", "h2", "e1", "e2", "b", "m1", "m2")]))
    rownames(coefs) <- c("h1", "h2", "e1", "e2", "b", "m1", "m2")
    colnames(coefs) <- c("True", "Est.")
    coefs[3:4,] <- exp(coefs[3:4,])
    coefs
  }, include.rownames = TRUE, digits = 3)

  ## Monotherapy plots
  output$marginals <- renderPlot({
    plot(fit(), logScale = input$logScale)
  })

  ## Table with occupancy values and constructed response
  output$occuptable <- DT::renderDataTable({

    d1e1 <- (compE()$occupancy$d1 / exp(pars()["e1"]))
    d2e2 <- (compE()$occupancy$d2 / exp(pars()["e2"]))
    occp1 <- (1 / compE()$occupancy$occupancy - 1)^(1/pars()["h1"])
    occp2 <- (1 / compE()$occupancy$occupancy - 1)^(1/pars()["h2"])
    weight1 <- d1e1 * occp1
    weight2 <- d2e2 * occp2

    contrib1 <- compE()$occupancy$occupancy * weight1 * (pars()["m1"] - pars()["b"])
    contrib2 <- compE()$occupancy$occupancy * weight2 * (pars()["m2"] - pars()["b"])
    baseline <- pars()["b"]

    ## Only valid if Hill coefficients are equal
    ## trueOcc <- with(comp()$occupancy, (d1/exp(pars()["e1"]) + d2/exp(pars()["e2"]))^pars()["h1"] /
    ##                                   (1 + (d1/exp(pars()["e1"]) + d2/exp(pars()["e2"]))^pars()["h1"]))

    printTable <- data.frame(comp()$occupancy[, c("d1", "d2")],
                             "Occupancy" = compE()$occupancy$occupancy,
                             "Weight1" = weight1,
                             "Weight2" = weight2,
                             "Baseline" = rep(baseline, nrow(compE()$occupancy)),
                             "Contrib1" = contrib1,
                             "Contrib2" = contrib2,
                             "Response" = comp()$response[-1])
    printTable <- printTable[order(abs(printTable$d1), decreasing = TRUE),]

    dat <- datatable(printTable,
                     options = list(pageLength = nrow(printTable),
                                    searching = FALSE),
                     rownames = FALSE,
                     selection = "single") %>%
      formatStyle(c("d1", "Weight1", "Contrib1"), color = "blue") %>%
      formatStyle(c("d2", "Weight2", "Contrib2"), color = "green") %>%
      formatRound(c("Occupancy", "Response", "Weight1", "Weight2",
                    "Contrib1", "Contrib2", "Baseline"), digits = 4) %>%
      formatRound(c("d1", "d2"), digits = 5)

  })


  output$isobologram <- renderPlot({

    surfaceFit <- list("data" = dr(), fitResult = fit(),
                       "null_model" = if (input$null == "stdloewe") "loewe" else input$null)
    class(surfaceFit) <- "ResponseSurface"
    isobologram(surfaceFit, logScale = input$logScale)

  })

  ## Response surface plot
  output$surface <-
    renderRglwidget({
      plotResponseSurface(data = dr(), fitResult = fit(), logScale = input$logScale,
                          null_model = if (input$null == "stdloewe") "loewe" else input$null,
                          legend = FALSE, colorBy = compE()$occupancy,
                          breaks = c(0, 0.25, 0.5, 0.75, 1),
                          plotfun = median,
                          colorPalette = c("#EFF3FF", "#BDD7E7", "#6BAED6", "#2171B5"))
      scene1 <- scene3d()
      close3d()
      rglwidget(scene1)
    })

})
