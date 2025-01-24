import os
import unittest
from unittest.mock import MagicMock, patch
import numpy as np
from melektrodica.collector import DataParameters, initialize


class TestDataParametersMarkdown(unittest.TestCase):
    def setUp(self):
        # Configurar un mock para 'writer' (simula el registrador de mensajes)
        self.writer = MagicMock()
        self.test_file = "test_parameters.md"

        # Crear un directorio de prueba, si no existe
        self.test_directory = "test_data"
        os.makedirs(self.test_directory, exist_ok=True)

        # Crear un archivo de prueba basado en los datos proporcionados
        self.test_file_path = os.path.join(self.test_directory, self.test_file)
        with open(self.test_file_path, "w") as f:
            f.write("""\
                        | Parameters                      | Variables                    | Value | Units |
                        |:--------------------------------|:-----------------------------|------:|:-----:|
                        | Electrode                       | Anode                        |  True |       |
                        | Physics parameters              | Temperature                  |    23 |   C   |
                        |                                 | Initial potential            |   0.0 |   V   |
                        |                                 | Final potential              |   0.5 |   V   |
                        |                                 | Step potential               | 0.001 |   V   |
                        | Continuous Stirred-Tank Reactor | Concentration = f(potential) | False |       |
                        |                                 | Catalyst Active surface area |   1.0 |  cm2  |
                        |                                 | Volumetric flux              |   1.0 |  L/s  |
                        | Rate Constants                  | --------------------------   |  ---- | ----  |
                        | Pre-exponential                 | A                            |     1 |       |
                        |                                 | j*                           |  True |       |
                        |                                 | j* (value)                   |   500 | A/cm2 |
                        | Transition state theory         | kappa * k_B * T^m / h        | False |       |
                        |                                 | kappa                        |     1 |       |
                        |                                 | m                            |     1 |       |
                        | Experimental                    | k_f, k_b                     | False |       |
                        | Thermochemical                  | Ga, DG_reaction              |  True |       |
                        |                                 | DG_reaction                  | False |       |
                        |                                 | G_formation                  |  True |       |
                        """)

    def tearDown(self):
        # Eliminar archivos y directorios temporales creados durante las pruebas
        if os.path.exists(self.test_file_path):
            os.remove(self.test_file_path)
        if os.path.exists(self.test_directory):
            os.rmdir(self.test_directory)

    @patch("melektrodica.Tool.unit_conversion", return_value=296.15)  # Conversión 23°C -> Kelvin
    @patch("melektrodica.Collector.raw_data")
    def test_initialization(self, mock_raw_data, mock_unit_conversion):
        # Simular los resultados de Collector.raw_data
        mock_header = ["Parameters", "Variables", "Value", "Units"]
        mock_data = np.array([
            ["Electrode", "Anode", "True", ""],
            ["Physics parameters", "Temperature", "23", "C"],
            ["", "Initial potential", "0.0", "V"],
            ["", "Final potential", "0.5", "V"],
            ["", "Step potential", "0.001", "V"],
            ["Continuous Stirred-Tank Reactor", "Concentration = f(potential)", "False", ""],
            ["", "Catalyst Active surface area", "1.0", "cm2"],
            ["", "Volumetric flux", "1.0", "L/s"],
            ["Rate Constants", "A", "1", ""],
            ["", "j*", "True", ""],
            ["", "j* (value)", "500", "A/cm2"],
            ["Transition state theory", "kappa * k_B * T^m / h", "False", ""],
            ["", "kappa", "1", ""],
            ["", "m", "1", ""],
            ["Experimental", "k_f, k_b", "False", ""],
            ["Thermochemical", "Ga, DG_reaction", "True", ""],
            ["", "DG_reaction", "False", ""],
            ["", "G_formation", "True", ""],
        ])
        mock_raw_data.return_value = (mock_header, mock_data)

        # Inicializar el objeto y verificar que funciona correctamente
        params = DataParameters(self.test_file_path, self.writer)

        # Verificar las propiedades inicializadas
        self.assertAlmostEqual(params.temperature, 296.15, places=2)  # Temperatura en Kelvin
        self.assertEqual(
            params.potential.tolist(),
            list(np.arange(0.0, 0.5 + 0.001, 0.001))  # Potenciales esperados [0.0, 0.001, ..., 0.5]
        )
        self.assertEqual(params.anode, True)
        self.assertFalse(params.experimental)
        self.assertTrue(params.thermochemical)
        self.assertEqual(params.pre_exponential, 1.0)
        self.assertEqual(params.js_value, 500.0)
        self.assertEqual(params.cstr, False)

        # Verificar mensajes enviados al writer mock
        self.writer.message.assert_any_call(f"Reading parameters from file: {self.test_file_path}")
        self.writer.message.assert_any_call("Thermochemical part details:\n\tDG_reaction: False\n\tG_formation: True")

        # Verificar mocks
        mock_unit_conversion.assert_called_once_with("Temperature", 23, "C", "K")
        mock_raw_data.assert_called_once_with(self.test_file_path, writer=self.writer)

    def test_initialize(self):
        # Verificar que la función 'initialize' funciona correctamente
        values = ["True", "296.15", "0.0", "0.5", "0.001", "True", "1.0", "1.0", "1", "True", "500", "False", "1", "1",
                  "False", "True", "False", "True"]
        lista = ["Anode", "Temperature", "Initial potential", "Final potential", "Step potential",
                 "Continuous Stirred-Tank Reactor", "Catalyst Active surface area", "Volumetric flux",
                 "A", "j*", "j* (value)", "kappa * k_B * T^m / h", "kappa", "m", "k_f, k_b",
                 "Ga, DG_reaction", "DG_reaction", "G_formation"]

        self.assertTrue(initialize(values, lista, "Anode"))
        self.assertTrue(initialize(values, lista, "j*"))
        self.assertFalse(initialize(values, lista, "Nonexistent parameter"))
        self.assertTrue(initialize(values, lista, "Continuous Stirred-Tank Reactor"))



if __name__ == "__main__":
    unittest.main()
