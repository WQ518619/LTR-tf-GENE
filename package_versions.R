# package_versions.R
library(here)

# 加载 step0_packages.R 并运行加载函数，获取包列表
source(here("R", "step0_packages.R"))
required_pkgs <- load_required_packages()  # 返回包名向量，但因为 invisible 不会打印

# 提取已加载包的版本信息
session_pkgs <- sessioninfo::session_info()$packages

# 只保留我们需要的包
pkg_info <- session_pkgs[session_pkgs$package %in% required_pkgs,
                         c("package", "loadedversion")]

# 保存结果
dir.create(here("results"), showWarnings = FALSE)
write.csv(pkg_info, here("results", "package_versions.csv"), row.names = FALSE)

message("Package versions saved to results/package_versions.csv")