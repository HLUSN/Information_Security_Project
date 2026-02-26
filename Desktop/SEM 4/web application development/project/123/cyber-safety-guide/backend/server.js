const express = require('express');
const cors = require('cors');
const mongoose = require('mongoose');

const app = express();

// Configure CORS to allow your Vercel frontend
const corsOptions = {
  origin: ['https://23-pi-rouge.vercel.app', 'http://localhost:3000'],
  methods: ['GET', 'POST', 'PUT', 'DELETE', 'OPTIONS'],
  allowedHeaders: ['Content-Type', 'Authorization'],
  credentials: true
};

app.use(cors(corsOptions));
app.use(express.json());

// Connect to MongoDB
mongoose.connect('mongodb+srv://admin:admin123@cluster0.necl9ae.mongodb.net/?retryWrites=true&w=majority&appName=Cluster0');

// Define Guide schema and model
const guideSchema = new mongoose.Schema({
  title: String,
  content: String,
  gifUrl: String,
});

const Guide = mongoose.model('Guide', guideSchema);

// Get all guides
app.get('/api/guides', async (req, res) => {
  try {
    const guides = await Guide.find();
    res.json(guides.map(g => ({ ...g.toObject(), id: g._id })));
  } catch (err) {
    res.status(500).json({ error: err.message });
  }
});

// Add a new guide
app.post('/api/guides', async (req, res) => {
  try {
    const guide = new Guide(req.body);
    const saved = await guide.save();
    res.json({ ...saved.toObject(), id: saved._id });
  } catch (err) {
    res.status(500).json({ error: err.message });
  }
});

// Update a guide
app.put('/api/guides/:id', async (req, res) => {
  try {
    const updated = await Guide.findByIdAndUpdate(req.params.id, req.body, { new: true });
    res.json({ ...updated.toObject(), id: updated._id });
  } catch (err) {
    res.status(500).json({ error: err.message });
  }
});

// Delete a guide
app.delete('/api/guides/:id', async (req, res) => {
  try {
    await Guide.findByIdAndDelete(req.params.id);
    res.json({ success: true });
  } catch (err) {
    res.status(500).json({ error: err.message });
  }
});

app.listen(8080, () => console.log('Server running on http://localhost:8080'));