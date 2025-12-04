import { cn } from "@/lib/utils";
import { Outlet, useLocation, useNavigate } from "react-router";
import { Tabs, TabsList, TabsTrigger } from "./ui/tabs";
import { Toaster } from "./ui/sonner";
import { Text } from "./ui/text";
import { GithubIcon } from "lucide-react";

export const Layout = () => {
  const navigate = useNavigate();
  const location = useLocation();

  const currentPath = location.pathname.substr(1);

  const handleChange = (val: string) => {
    if (val === "") navigate("/");
    if (val === "proposal") navigate("/proposal");
    if (val === "midterm") navigate("/midterm");
    if (val === "report") navigate("/report");
  };

  return (
    <div
      className={cn(
        "absolute w-screen min-h-screen top-0 left-0 overflow-y-auto m-0 p-0 scroll-smooth",
        "bg-neutral-900 pt-20",
        "flex flex-col items-center",
      )}
    >
      <div className="print-hide z-100 w-full fixed top-0 pt-4 pb-8 flex justify-center bg-linear-to-b from-neutral-900 to-neutral-900/0">
        <Tabs value={currentPath} onValueChange={handleChange}>
          <TabsList>
            <TabsTrigger value="">Home</TabsTrigger>
            <TabsTrigger value="proposal">Proposal</TabsTrigger>
            <TabsTrigger value="midterm">Midterm Checkpoint</TabsTrigger>
            <TabsTrigger value="report">Final Report</TabsTrigger>
          </TabsList>
        </Tabs>
        <a href="https://github.com/i3ta/brain-storm" target="_blank">
          <GithubIcon className="absolute top-6 right-8 cursor-pointer hover:opacity-50 transition-all" />
        </a>
      </div>
      <Outlet />
      <div className="print-padding" />
      <div
        className={cn(
          "print-hide",
          "w-11/12 max-w-5xl p-4 mt-16 mb-8 rounded-xl border bg-neutral-800/50",
          "flex flex-col items-center",
        )}
      >
        <Text size="c">
          Written By Team 24 for Georgia Tech CS 4641, Fall 2025
        </Text>
        <Text size="c">
          Copyright 2025 © Ritika Calyanakoti, Aaron Fan, Hannah Hsiao, Aaron
          Hung, Sagarika Srinivasan
        </Text>
      </div>
      <Toaster richColors />
    </div>
  );
};
