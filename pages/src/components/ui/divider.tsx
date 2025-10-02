import { cn } from "@/lib/utils";
import type { ComponentProps } from "react";

export const Divider = ({ className, ...props }: ComponentProps<"div">) => {
  return (
    <div
      className={cn("w-full mb-4 border-b border-neutral-500", className)}
      {...props}
    />
  );
};
