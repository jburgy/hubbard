const express = require('express');
const morgan = require('morgan');
const app = express();
const port = 3000;

express.static.mime.types['wasm'] = 'application/wasm';

app.use(morgan('common'));
app.use(express.static(process.cwd()));

app.listen(port, () => console.log(`server listening on port ${port}`));