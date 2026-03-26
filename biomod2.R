library(biomod2)
# ==========================================
# 1. 准备物种存在点 (XY坐标)
# ==========================================
# 从你清洗好的数据中，提取出这个物种的纯坐标
sp_name <- "Abaeis_nicippe"
myRespName <- sp_name

myRespXY <- pixel_level_df %>%
  filter(species == sp_name) %>%
  dplyr::select(x, y) # 只保留经纬度坐标

# 创建一个全是 1 的向量，代表这些坐标点都是“存在(Presence)”的
myResp <- rep(1, nrow(myRespXY)) 


# ==========================================
# 2. 准备环境变量栅格 (也就是我们之前的多层图层)
# ==========================================
# 记得把 landuse 声明为分类变量！这对 biomod2 至关重要
# 假设 env_cont 和 landuse_ras 已经在你的环境中
myExpl <- c(env_cont, landuse_ras)

# 告诉 R 这是一个分类变量 (Factor)
myExpl$landuse <- as.factor(myExpl$landuse)


# ==========================================
# 3. 让 biomod2 格式化数据 (重头戏！)
# ==========================================
myBiomodData <- BIOMOD_FormatingData(
  resp.var = myResp,           # 物种存在状态 (全是 1)
  expl.var = myExpl,           # 环境变量栅格图层
  resp.xy = myRespXY,          # 物种所在的 XY 坐标
  resp.name = myRespName,      # 物种名字
  PA.nb.rep = 1,               # 生成几组伪缺席点 (测试先设1组)
  PA.nb.absences = 1000,       # 随机生成 1000 个伪缺席点 (背景点)
  PA.strategy = 'random'       # 在图层范围内随机撒点
)

# 打印出来看看 biomod2 是否成功接管了数据！
myBiomodData
plot(myBiomodData) # 还能直接画出点位图！