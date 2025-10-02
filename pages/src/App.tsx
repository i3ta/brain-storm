import { BrowserRouter, Route, Routes } from "react-router";
import { Home } from "./pages/homepage/homepage";
import { Proposal } from "./pages/proposal/proposal";
import { Layout } from "./components/layout";

function App() {
  return (
    <BrowserRouter>
      <Routes>
        <Route element={<Layout />}>
          <Route path="/" element={<Home />} />
          <Route path="/proposal" element={<Proposal />} />
        </Route>
      </Routes>
    </BrowserRouter>
  );
}

export default App;
