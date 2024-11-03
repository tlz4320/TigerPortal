anno_barplot2 <- function (x, baseline = 0, which = c("column", "row"), border = TRUE, 
                           bar_width = 0.6, beside = FALSE, attach = FALSE, gp = gpar(fill = "#CCCCCC"), 
                           ylim = NULL, extend = 0.05, axis = TRUE, axis_param = default_axis_param(which), 
                           add_numbers = FALSE, numbers_gp = gpar(fontsize = 8), numbers_rot = ifelse(which == 
                                                                                                        "column", 45, 0), numbers_offset = unit(2, "mm"), width = NULL, 
                           height = NULL, ...) 
{
  other_args = list(...)
  if (length(other_args)) {
    if ("axis_gp" %in% names(other_args)) {
      stop_wrap("`axis_gp` is removed from the arguments. Use `axis_param = list(gp = ...)` instead.")
    }
    if ("axis_side" %in% names(other_args)) {
      stop_wrap("`axis_side` is removed from the arguments. Use `axis_param = list(side = ...)` instead.")
    }
    if ("axis_direction" %in% names(other_args)) {
      stop_wrap("`axis_direction` is not supported any more.")
    }
  }
  if (inherits(x, "list")) 
    x = do.call("cbind", x)
  if (inherits(x, "data.frame")) 
    x = as.matrix(x)
  if (inherits(x, "matrix")) {
    sg = apply(x, 1, function(xx) all(sign(xx) %in% c(1, 
                                                      0)) || all(sign(xx) %in% c(-1, 0)))
    if (!all(sg)) {
      stop_wrap("Since `x` is a matrix, the sign of each row should be either all positive or all negative.")
    }
  }
  labels_format = attr(x, "labels_format")
  if (is.null(dim(x))) 
    x = matrix(x, ncol = 1)
  nc = ncol(x)
  if (missing(gp)) {
    gp = gpar(fill = grey(seq(0, 1, length.out = nc + 2))[-c(1, 
                                                             nc + 2)])
  }
  if (beside) {
    data_scale = range(x, na.rm = TRUE)
  }
  else {
    data_scale = range(rowSums(x, na.rm = TRUE), na.rm = TRUE)
  }
  if (data_scale[1] == data_scale[2]) 
    data_scale[2] = data_scale[1] + .Machine$double.eps * 
    1.1
  if (!is.null(ylim)) 
    data_scale = ylim
  if (baseline == "min") {
    data_scale = data_scale + c(0, extend) * (data_scale[2] - 
                                                data_scale[1])
    baseline = min(x, na.rm = TRUE)
  }
  else if (baseline == "max") {
    data_scale = data_scale + c(-extend, 0) * (data_scale[2] - 
                                                 data_scale[1])
    baseline = max(x, na.rm = TRUE)
  }
  else {
    if (is.numeric(baseline)) {
      if (baseline == 0 && all(abs(rowSums(x, na.rm = TRUE) - 
                                   1) < 1e-06) && !beside) {
        data_scale = c(0, 1)
      }
      else if (baseline <= data_scale[1]) {
        data_scale = c(baseline, extend * (data_scale[2] - 
                                             baseline) + data_scale[2])
      }
      else if (baseline >= data_scale[2]) {
        data_scale = c(-extend * (baseline - data_scale[1]) + 
                         data_scale[1], baseline)
      }
      else {
        data_scale = data_scale + c(-extend, extend) * 
          (data_scale[2] - data_scale[1])
      }
    }
  }
  ef = function() NULL
  if (is.null(ComplexHeatmap:::.ENV$current_annotation_which)) {
    which = match.arg(which)[1]
    dev.null()
    ef = dev.off2
  }
  else {
    which = ComplexHeatmap:::.ENV$current_annotation_which
  }
  on.exit(ef())
  anno_size = ComplexHeatmap:::anno_width_and_height(which, width, height, 
                                                     unit(1, "cm"))
  if (nc == 1) {
    gp = ComplexHeatmap:::recycle_gp(gp, nrow(x))
  }
  else {
    gp = ComplexHeatmap:::recycle_gp(gp, nc)
  }
  value = x
  attr(value, "labels_format") = labels_format
  if (ncol(value) == 1) {
    if (add_numbers) {
      if (which == "column") {
        if (numbers_rot == 0) {
          extend = convertHeight(max_text_height(value, 
                                                 gp = numbers_gp) + numbers_offset + unit(2, 
                                                                                          "mm"), "mm", valueOnly = TRUE)/convertHeight(anno_size$height, 
                                                                                                                                       "mm", valueOnly = TRUE) * (data_scale[2] - 
                                                                                                                                                                    data_scale[1])
        }
        else {
          extend = convertHeight(sin(numbers_rot/180 * 
                                       pi) * max_text_width(value, gp = numbers_gp) + 
                                   numbers_offset + unit(4, "mm"), "mm", valueOnly = TRUE)/convertHeight(anno_size$height, 
                                                                                                         "mm", valueOnly = TRUE) * (data_scale[2] - 
                                                                                                                                      data_scale[1])
        }
        data_scale[2] = data_scale[2] + extend
      }
      else if (which == "row") {
        extend = convertWidth(cos(numbers_rot/180 * 
                                    pi) * max_text_width(value, gp = numbers_gp) + 
                                numbers_offset + unit(4, "mm"), "mm", valueOnly = TRUE)/convertWidth(anno_size$width, 
                                                                                                     "mm", valueOnly = TRUE) * (data_scale[2] - 
                                                                                                                                  data_scale[1])
        data_scale[2] = data_scale[2] + extend
      }
    }
  }
  axis_param = ComplexHeatmap:::validate_axis_param(axis_param, which)
  axis_grob = if (axis) 
    ComplexHeatmap:::construct_axis_grob(axis_param, which, data_scale, format = labels_format)
  else NULL
  row_fun = function(index, k = 1, N = 1) {
    n = length(index)
    if (axis_param$direction == "reverse") {
      value_origin = value
      value = data_scale[2] - value + data_scale[1]
      baseline = data_scale[2] - baseline + data_scale[1]
    }
    pushViewport(viewport(xscale = data_scale, yscale = c(0.5, 
                                                          n + 0.5)))
    if (ncol(value) == 1) {
      width = value[index] - baseline
      x_coor = width/2 + baseline
      grid.rect(x = x_coor, y = n - seq_along(index) + 
                  1, width = abs(width), height = 1 * bar_width, 
                default.units = "native", gp = subset_gp(gp, 
                                                         index))
      if (add_numbers) {
        if (axis_param$direction == "normal") {
          txt = value[index]
          if (!is.null(attr(value, "labels_format"))) {
            txt = attr(value, "labels_format")(value[index])
          }
          grid.text(txt, x = unit(baseline + width, 
                                  "native") + numbers_offset, y = n - seq_along(index) + 
                      1, default.units = "native", gp = subset_gp(numbers_gp, 
                                                                  index), just = c("left"), rot = numbers_rot)
        }
        else {
          txt = value_origin[index]
          if (!is.null(attr(value, "labels_format"))) {
            txt = attr(value, "labels_format")(value_origin[index])
          }
          grid.text(txt, x = unit(baseline + width, 
                                  "native") - numbers_offset, y = n - seq_along(index) + 
                      1, default.units = "native", gp = subset_gp(numbers_gp, 
                                                                  index), just = c("right"), rot = numbers_rot)
        }
      }
    }
    else {
      if (beside) {
        nbar = ncol(value)
        nr = length(index)
        for (i in seq_along(index)) {
          for (j in 1:nbar) {
            if (attach) {
              if (axis_param$direction == "normal") {
                grid.rect(x = baseline, y = nr - i + 
                            0.5 + (1 - bar_width)/2 + (nbar - 
                                                         j + 0.5)/nbar * bar_width, width = value[index[i], 
                                                                                                  j], height = 1/nbar * bar_width, just = c("left"), 
                          default.units = "native", gp = subset_gp(gp, 
                                                                   j))
              }
              else {
                grid.rect(x = baseline, y = nr - i + 
                            0.5 + (1 - bar_width)/2 + (nbar - 
                                                         j + 0.5)/nbar * bar_width, width = value[index[i], 
                                                                                                  j], height = 1/nbar * bar_width, just = c("right"), 
                          default.units = "native", gp = subset_gp(gp, 
                                                                   j))
              }
            }
            else {
              if (axis_param$direction == "normal") {
                grid.rect(x = baseline, y = nr - i + 
                            0.5 + (nbar - j + 0.5)/nbar, width = value[index[i], 
                                                                       j], height = 1/nbar * bar_width, just = c("left"), 
                          default.units = "native", gp = subset_gp(gp, 
                                                                   j))
              }
              else {
                grid.rect(x = baseline, y = nr - i + 
                            0.5 + (nbar - j + 0.5)/nbar, width = value[index[i], 
                                                                       j], height = 1/nbar * bar_width, just = c("right"), 
                          default.units = "native", gp = subset_gp(gp, 
                                                                   j))
              }
            }
          }
        }
      }
      else {
        for(index2 in index){
          value_tmp <- unlist(value[index2,])
          order_tmp <- order(value_tmp, decreasing = T)
          value_tmp <- sort(value_tmp, decreasing = T)
          for (i in seq_len(ncol(value))) {
            if (axis_param$direction == "normal") {
              width = abs(value_tmp[i])
              x_coor = sum(value_tmp[seq_len(i - 1), drop = FALSE]) + width/2
              grid.rect(x = x_coor, y = n - index2 + 
                          1, width = abs(width), height = 1 * bar_width, 
                        default.units = "native", gp = subset_gp(gp, order_tmp[i]))
            }
            else {
              width = value_origin[index2, i]
              x_coor = rowSums(value_origin[index2, seq_len(i - 
                                                              1), drop = FALSE]) + width/2
              x_coor = data_scale[2] - x_coor + data_scale[1]
              grid.rect(x = x_coor, y = n - seq_along(index2) + 
                          1, width = abs(width), height = 1 * bar_width, 
                        default.units = "native", gp = subset_gp(gp, 
                                                                 i))
            }
          }
        }
      }
    }
    if (axis_param$side == "top") {
      if (k > 1) 
        axis = FALSE
    }
    else if (axis_param$side == "bottom") {
      if (k < N) 
        axis = FALSE
    }
    if (axis) 
      grid.draw(axis_grob)
    if (border) 
      grid.rect(gp = gpar(fill = "transparent"))
    popViewport()
  }
  column_fun = function(index, k = 1, N = 1) {
    n = length(index)
    if (axis_param$direction == "reverse") {
      value_origin = value
      value = data_scale[2] - value + data_scale[1]
      baseline = data_scale[2] - baseline + data_scale[1]
    }
    pushViewport(viewport(yscale = data_scale, xscale = c(0.5, 
                                                          n + 0.5)))
    if (ncol(value) == 1) {
      height = value[index] - baseline
      y_coor = height/2 + baseline
      grid.rect(y = y_coor, x = seq_along(index), height = abs(height), 
                width = 1 * bar_width, default.units = "native", 
                gp = subset_gp(gp, index))
      if (add_numbers) {
        txt = value[index]
        if (!is.null(attr(value, "labels_format"))) {
          txt = attr(value, "labels_format")(value[index])
        }
        numbers_rot = numbers_rot%%360
        if (numbers_rot == 0) {
          grid.text(txt, x = seq_along(index), y = unit(baseline + 
                                                          height, "native") + numbers_offset, default.units = "native", 
                    gp = subset_gp(numbers_gp, index), just = c("bottom"))
        }
        else {
          grid.text(txt, x = seq_along(index), y = unit(baseline + 
                                                          height, "native") + numbers_offset, default.units = "native", 
                    gp = subset_gp(numbers_gp, index), just = c("left"), 
                    rot = numbers_rot)
        }
      }
    }
    else {
      if (beside) {
        nbar = ncol(value)
        nr = length(index)
        for (i in seq_along(index)) {
          for (j in 1:nbar) {
            if (attach) {
              if (axis_param$direction == "normal") {
                grid.rect(y = baseline, x = i - 0.5 + 
                            (1 - bar_width)/2 + (j - 0.5)/nbar * 
                            bar_width, height = value[index[i], 
                                                      j], width = 1/nbar * bar_width, just = c("bottom"), 
                          default.units = "native", gp = subset_gp(gp, 
                                                                   j))
              }
              else {
                grid.rect(y = baseline, x = i - 0.5 + 
                            (1 - bar_width)/2 + (j - 0.5)/nbar * 
                            bar_width, height = value[index[i], 
                                                      j], width = 1/nbar * bar_width, just = c("top"), 
                          default.units = "native", gp = subset_gp(gp, 
                                                                   j))
              }
            }
            else {
              if (axis_param$direction == "normal") {
                grid.rect(y = baseline, x = i - 0.5 + 
                            (j - 0.5)/nbar, height = value[index[i], 
                                                           j], width = 1/nbar * bar_width, just = c("bottom"), 
                          default.units = "native", gp = subset_gp(gp, 
                                                                   j))
              }
              else {
                grid.rect(y = baseline, x = i - 0.5 + 
                            (j - 0.5)/nbar, height = value[index[i], 
                                                           j], width = 1/nbar * bar_width, just = c("top"), 
                          default.units = "native", gp = subset_gp(gp, 
                                                                   j))
              }
            }
          }
        }
      }
      else {
        for (i in seq_len(ncol(value))) {
          if (axis_param$direction == "normal") {
            height = value[index, i]
            y_coor = rowSums(value[index, seq_len(i - 
                                                    1), drop = FALSE]) + height/2
            grid.rect(y = y_coor, x = seq_along(index), 
                      height = abs(height), width = 1 * bar_width, 
                      default.units = "native", gp = subset_gp(gp, 
                                                               i))
          }
          else {
            height = value_origin[index, i]
            y_coor = rowSums(value_origin[index, seq_len(i - 
                                                           1), drop = FALSE]) + height/2
            y_coor = data_scale[2] - y_coor + data_scale[1]
            grid.rect(y = y_coor, x = seq_along(index), 
                      height = abs(height), width = 1 * bar_width, 
                      default.units = "native", gp = subset_gp(gp, 
                                                               i))
          }
        }
      }
    }
    if (axis_param$side == "left") {
      if (k > 1) 
        axis = FALSE
    }
    else if (axis_param$side == "right") {
      if (k < N) 
        axis = FALSE
    }
    if (axis) 
      grid.draw(axis_grob)
    if (border) 
      grid.rect(gp = gpar(fill = "transparent"))
    popViewport()
  }
  if (which == "row") {
    fun = row_fun
  }
  else if (which == "column") {
    fun = column_fun
  }
  n = nrow(value)
  anno = AnnotationFunction(fun = fun, fun_name = "anno_barplot", 
                            which = which, width = anno_size$width, height = anno_size$height, 
                            n = n, data_scale = data_scale, var_import = list(value, 
                                                                              gp, border, bar_width, baseline, beside, attach, 
                                                                              axis, axis_param, axis_grob, data_scale, add_numbers, 
                                                                              numbers_gp, numbers_offset, numbers_rot))
  anno@subset_rule$value = subset_matrix_by_row
  if (ncol(value) == 1) {
    anno@subset_rule$gp = subset_gp
  }
  anno@subsettable = TRUE
  anno@extended = ComplexHeatmap:::update_anno_extend(anno, axis_grob, axis_param)
  return(anno)
}