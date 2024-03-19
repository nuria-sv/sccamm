library(latex2exp)

XX2 = X
YY2 = Y
# temporal dependence check
j = 80
X_id = XX2[XX2$id == j, ]
Y_id = YY2[YY2$id == j, ]

dim(X_id);dim(Y_id)

dev.new()
par(mfrow = c(1, 2))
image(cov(t(X_id[,-c(1,2)])), col = hcl.colors(12, "Purple-Blue", rev = TRUE), axes = FALSE)
# contour(cor(t(X_id[-c(1,2), ])), add = TRUE, drawlabels = FALSE)
axis(1, at = seq(0, 1, length=length(X_id$id)), labels = X_id$time)
axis(2, at = seq(0, 1, length=length(X_id$id)), labels = X_id$time)
title(main = paste0("measurement cov in X for i = ", j), font.main = 4)

image(cov(t(Y_id[-c(1,2), ])), col = hcl.colors(12, "Purple-Blue", rev = TRUE), axes = FALSE)
axis(1, at = seq(0, 1, length=length(Y_id$id)), labels = Y_id$time)
axis(2, at = seq(0, 1, length=length(Y_id$id)), labels = Y_id$time)
title(main = paste0("measurement cov in Y for i = ", j), font.main = 4)


# plots after cca
z_x = predict(res_toscca$me_x, newdata = data.frame(time = X_id$time, id = X_id$id), allow.new.levels = TRUE, re.form = NULL)
z_y = predict(res_toscca$me_y, newdata = data.frame(time = Y_id$time, id = Y_id$id), allow.new.levels = TRUE, re.form = NULL)

lv_x = data.frame(cbind(X_id[, c(1,2)], lv_x = as.matrix(z_x)))
lv_y = data.frame(cbind(Y_id[, c(1,2)], lv_y = as.matrix(z_y)))

# lv_x_id = lv_x[lv_x$id == j, ]
# lv_y_id = lv_y[lv_y$id == j, ]

new_x = data.frame(cbind(X_id[, c(1,2)], lv_x[,-c(1,2)]%*%t((res_toscca$a))))
new_y = data.frame(cbind(Y_id[, c(1,2)], lv_y[,-c(1,2)]%*%t((res_toscca$b))))

# new_x_id = new_x[new_x$id == j, ]
# new_y_id = new_y[new_y$id == j, ]

dev.new()
par(mfrow = c(1, 2))
image(cor(t(new_x[,-c(1,2)])), col = hcl.colors(12, "Purple-Blue", rev = TRUE), axes = FALSE)
# contour(cor(t(X_id[-c(1,2), ])), add = TRUE, drawlabels = FALSE)
axis(1, at = seq(0, 1, length=length(X_id$id)), labels = X_id$time)
axis(2, at = seq(0, 1, length=length(X_id$id)), labels = X_id$time)
title(main = paste0("measurement cov in ", expression(hat(X)), " for i = ", j), font.main = 4)

image(cov(t(new_y[,-c(1,2)])), col = hcl.colors(12, "Purple-Blue", rev = TRUE), axes = FALSE)
axis(1, at = seq(0, 1, length=length(Y_id$id)), labels = Y_id$time)
axis(2, at = seq(0, 1, length=length(Y_id$id)), labels = Y_id$time)
title(main = paste0("measurement cov in ", expression(hat(Y)), " for i = ", j), font.main = 4)

error = X_id[,-c(1,2)] - new_x[,-c(1,2)]
image(cor(t(error)), col = hcl.colors(12, "Purple-Blue", rev = TRUE), axes = FALSE)
axis(1, at = seq(0, 1, length=length(Y_id$id)), labels = Y_id$time)
axis(2, at = seq(0, 1, length=length(Y_id$id)), labels = Y_id$time)
title(main = paste0("measurement cov in residuals for i = ", j), font.main = 4)

z_x = predict(res_toscca$me_x, newdata = data.frame(time = XX2$time, id = XX2$id), allow.new.levels = TRUE, re.form = NULL)
z_y = predict(res_toscca$me_y, newdata = data.frame(time = YY2$time, id = YY2$id), allow.new.levels = TRUE, re.form = NULL)
