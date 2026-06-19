import struct
import array
import os

def convert_txt_to_binary(input_filename, output_filename):
    """Liest eine Textdatei mit Floats und speichert sie im spezifischen Binärformat."""
    
    if not os.path.exists(input_filename):
        print(f"Fehler: Die Datei '{input_filename}' wurde nicht gefunden.")
        return

    print(f"Lese Daten aus '{input_filename}'...")
    
    # 1. Daten einlesen
    data = []
    with open(input_filename, "r", encoding="utf-8") as infile:
        for line in infile:
            line = line.strip()
            if line:
                data.append(float(line))

    # 2. Metadaten für den Header definieren
    bytes_of_data_type = 4          # 4 Bytes für 32-bit Float (Standard für solche Formate)
    size_of_vectors = 1             # Annahme: Alles ist eine einzige, lange Zeitreihe
    dimension = len(data)           # Die Dimension ist die gesamte Anzahl an Zeilen/Datenpunkten
    
    print(f"Schreibe Binärdatei '{output_filename}'...")
    print(f" - Bytes pro Wert: {bytes_of_data_type}")
    print(f" - Anzahl Vektoren: {size_of_vectors}")
    print(f" - Dimension: {dimension}")

    # 3. Binärdatei schreiben
    with open(output_filename, "wb") as outfile:
        # Header schreiben (jeweils als 4-Byte Integer, 'i' steht für C-Type int)
        # '<' erzwingt Little-Endian (Standard bei den meisten Systemen und C/C++ Programmen)
        outfile.write(struct.pack('<i', bytes_of_data_type))
        outfile.write(struct.pack('<i', size_of_vectors))
        outfile.write(struct.pack('<i', dimension))

        # Daten schreiben
        # 'f' konvertiert die Python-Floats in 4-Byte (32-bit) C-Floats
        # Dies ist wesentlich schneller und speichereffizienter als struct.pack in einer Schleife
        float_array = array.array('f', data)
        float_array.tofile(outfile)

    print("\nFertig! Konvertierung erfolgreich abgeschlossen.")

# Skript ausführen
convert_txt_to_binary("test.txt", "test.data_new")