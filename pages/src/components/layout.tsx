import { cn } from "@/lib/utils";
import { Outlet, useLocation, useNavigate } from "react-router";
import { Tabs, TabsList, TabsTrigger } from "./ui/tabs";
import { Toaster } from "./ui/sonner";

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
      <div className="w-full fixed top-4 flex justify-center">
        <Tabs value={currentPath} onValueChange={handleChange}>
          <TabsList>
            <TabsTrigger value="">Home</TabsTrigger>
            <TabsTrigger value="proposal">Proposal</TabsTrigger>
            <TabsTrigger value="midterm">Midterm Checkpoint</TabsTrigger>
            <TabsTrigger value="report">Final Report</TabsTrigger>
          </TabsList>
        </Tabs>
      </div>
      <Outlet />
      <Toaster richColors />
    </div>
  );
};
