(function () {
  // 年份
  var el = document.getElementById("year");
  if (el) {
    el.textContent = String(new Date().getFullYear());
  }

  // 生成闪烁星星
  function initStars() {
    const starsContainer = document.querySelector(".stars-bg");
    if (!starsContainer) return;

    const starCount = Math.floor(window.innerWidth / 50); // 根据屏幕宽度动态设置星星数量
    const maxStars = Math.min(starCount, 100); // 最多100个星星

    for (let i = 0; i < maxStars; i++) {
      const star = document.createElement("div");
      star.classList.add("star", "animate");

      // 随机位置
      const x = Math.random() * 100;
      const y = Math.random() * 100;
      star.style.left = x + "%";
      star.style.top = y + "%";

      // 随机大小（1-3px）
      const size = Math.random() * 2 + 1;
      star.style.width = size + "px";
      star.style.height = size + "px";

      // 随机闪烁时长 (2-4秒)
      const duration = Math.random() * 2 + 2;
      star.style.setProperty("--duration", duration + "s");

      // 随机延迟开始
      const delay = Math.random() * duration;
      star.style.animationDelay = delay + "s";

      // 循环闪烁时重新随机位置和大小
      star.addEventListener("animationiteration", function (event) {
        if (event.animationName !== "twinkle") return;
        const x2 = Math.random() * 100;
        const y2 = Math.random() * 100;
        star.style.left = x2 + "%";
        star.style.top = y2 + "%";
        const newSize = Math.random() * 2 + 1;
        star.style.width = newSize + "px";
        star.style.height = newSize + "px";
      });

      starsContainer.appendChild(star);
    }
  }

  // 初始化星星
  initStars();

  // 响应式重新生成星星
  let resizeTimer;
  window.addEventListener("resize", function () {
    clearTimeout(resizeTimer);
    resizeTimer = setTimeout(function () {
      document.querySelector(".stars-bg").innerHTML = "";
      initStars();
    }, 300);
  });
})();
