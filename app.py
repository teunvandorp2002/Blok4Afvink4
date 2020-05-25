from flask import Flask, render_template, request
import sequence_processor as seq

app = Flask(__name__)


@app.route('/', methods=['POST', 'GET'])
def input_():
    sequence = str()
    st = str()
    if request.method == 'POST':
        sequence = request.form["sequence"]
        st, sequence = seq.sequence_type(sequence)
    if st == 'DNA':
        dna_attributes = seq.dna_processor(sequence)
        return render_template("output.html", dna=dna_attributes)
    if st == 'Eiwit':
        seq.blast(sequence)
        blasts = seq.parse()
        return render_template("output.html", blasts=blasts)
    return render_template("Input.html", sequence=sequence, seq_type=st)


if __name__ == '__main__':
    app.run()
