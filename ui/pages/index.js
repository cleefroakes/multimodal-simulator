import React from 'react';

export default function Home() {
  return (
    <div>
      <h1>Multimodal Simulator</h1>
      <input type="text" placeholder="Enter text or upload image" />
      <button onClick={() => console.log('Simulate')}>Simulate</button>
    </div>
  );
}