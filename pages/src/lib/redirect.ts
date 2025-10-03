import { useEffect } from "react";
import { useNavigate, useLocation } from "react-router";

export const RedirectHandler = () => {
  const navigate = useNavigate();
  const location = useLocation();

  useEffect(() => {
    const hash = window.location.hash;
    const params = new URLSearchParams(window.location.search);
    const redirectPath = params.get("redirect");

    if (redirectPath) {
      console.log("redirect", redirectPath);
      navigate("/" + redirectPath, { replace: true });
    }

    if (hash) {
      console.log("going to hash", hash);
      const id = hash.replace("#", "");
      let attempts = 0;
      const maxAttempts = 20;

      const scrollToElement = () => {
        const el = document.getElementById(id);
        if (el) {
          el.scrollIntoView({ behavior: "smooth" });
          console.log("Successfully scrolled to", id);
        } else if (attempts < maxAttempts) {
          attempts++;
          setTimeout(scrollToElement, 100);
        } else {
          console.warn("Could not find element with id:", id);
        }
      };

      setTimeout(scrollToElement, 300);
    }
  }, [navigate, location.pathname, location.hash]);

  return null;
};
