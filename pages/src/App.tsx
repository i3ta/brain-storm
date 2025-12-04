import { BrowserRouter, Route, Routes } from "react-router";
import { Home } from "./pages/homepage/homepage";
import { Proposal } from "./pages/proposal/proposal";
import { Layout } from "./components/layout";
import { RedirectHandler } from "./lib/redirect";
import { Midterm } from "./pages/midterm/midterm";
import { Report } from "./pages/report/report";

const basename = "/brain-storm/";

function App() {
  return (
    <BrowserRouter basename={basename}>
      <RedirectHandler />
      <Routes>
        <Route element={<Layout />}>
          <Route path="/" element={<Home />} />
          <Route path="/proposal" element={<Proposal />} />
          <Route path="/midterm" element={<Midterm />} />
          <Route path="/report" element={<Report />} />
        </Route>
      </Routes>
    </BrowserRouter>
  );
}

export default App;
