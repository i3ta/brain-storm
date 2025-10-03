import { copy } from "@/lib/copy";
import { cn } from "@/lib/utils";
import type { ComponentProps } from "react";

export interface TextProps extends ComponentProps<"p"> {
  size?: "t1" | "t2" | "h1" | "h2" | "h3" | "p" | "c";
  id?: string;
}

export const Text = ({
  children,
  size = "p",
  className,
  id,
  ...props
}: TextProps) => {
  const sizeClasses: Record<NonNullable<TextProps["size"]>, string> = {
    t1: "text-[96px] font-bold text-transparent bg-clip-text bg-gradient-to-r from-[#11998e] to-[#38ef7d] leading-none pb-8",
    t2: "text-[64px] font-bold",
    h1: "text-[32px] font-bold pt-8",
    h2: "text-[28px] font-bold pt-8",
    h3: "text-[20px] font-bold",
    p: "text-[18px] font-normal pb-2",
    c: "text-[12px] font-normal pb-2",
  };

  const handleCopyLink = () => {
    if (!id) return;
    const fullUrl = window.location.href;
    copy(`${fullUrl}#${id}`);
  };

  return (
    <>
      <p
        className={cn(
          "relative text-white inline",
          sizeClasses[size],
          className,
        )}
        id={id}
        {...props}
      >
        {children}
        {id && size != "p" && (
          <span
            className={cn(
              "absolute top-0 -left-2.5 -translate-x-full transition-all cursor-pointer",
              sizeClasses[size],
              "text-neutral-700 hover:text-neutral-500",
              className,
            )}
            onClick={handleCopyLink}
          >
            #
          </span>
        )}
      </p>
    </>
  );
};
