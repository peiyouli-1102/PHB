# PHB

个人主页与博客（纯静态 HTML/CSS/JS）。https://peiyouli-1102.github.io/PHB/

## 本地预览

用浏览器直接打开仓库根目录下的 `index.html`，或使用任意静态文件服务（例如 `python -m http.server`）后访问 `http://localhost:8000/`。

## GitHub Pages 部署（依据 PLAN 阶段 5）

本站全部使用**相对路径**，适用于「用户主页」或「项目站」子路径；站点文件在仓库根目录。

### 1. 在 GitHub 新建仓库

在 GitHub 上创建一个空仓库（**不要**勾选「Add a README」可避免首次推送冲突）。记下 HTTPS 地址，例如：

`https://github.com/<用户名>/<仓库名>.git`

### 2. 本地推送

在本仓库根目录执行：

```bash
git remote add origin https://github.com/<用户名>/<仓库名>.git
git push -u origin main
```

若远程仓库初始化时自带 README，可先拉取再推送：

```bash
git pull origin main --allow-unrelated-histories
# 解决若有冲突后
git push -u origin main
```

### 3. 启用 Pages

打开仓库 **Settings → Pages**：

- **Source**：Deploy from a branch  
- **Branch**：`main`，文件夹 **`/ (root)`**（站点入口为根目录的 `index.html`）

保存后等待一两分钟，页面顶部会显示站点 URL（若为用户主页仓库 `<用户名>.github.io`，否则为 `https://<用户名>.github.io/<仓库名>/`）。

### 4. 上线后自检（PLAN 验收）

用无痕窗口打开线上首页，依次验证：**首页 → 博客列表 → 第一篇 → 下载附件 → 返回主页**，确认无 404。
