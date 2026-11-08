import numpy as np
import scipy.linalg as spl
import statistics
import functools
import matplotlib.pyplot as plt
import difflib
hadamard_matrix = spl.hadamard(16)
SF_ONE = [0, 1]
SF_ZERO = [1, 0]
SF_BIT_LEN = len(SF_ONE)

# Header CRC-8. 'legacy' (default): polynomial 0xD5, register initialised to 0, the 24 header bits fed one bit per
# byte through the byte table (a [32,24] code of minimum distance 2); every printed result and the SDR chip file use
# it. 'standard': IEEE 802.15.6-2012 10.6.1, G(x) = x^8 + x^7 + x^3 + x^2 + 1 (0x8D), register initialised to all
# ones, bits processed in transmission order, output complemented (minimum distance 3 on the 24 header bits).
HEADER_CRC_MODE = 'legacy'


def header_crc_bits(bits, mode=None):
    """8-character bit string of the header CRC over a list of 0/1 header bits in transmission order."""
    mode = HEADER_CRC_MODE if mode is None else mode
    if mode == 'legacy':
        return format(CRC8().compute([int(b) for b in bits]), '08b')
    crc = 0xFF
    for b in bits:
        crc ^= int(b) << 7
        crc = ((crc << 1) ^ 0x8D) & 0xFF if crc & 0x80 else (crc << 1) & 0xFF
    return format(crc ^ 0xFF, '08b')


class CRC8:

    def __init__(self, polynomial=0xD5):
        self.polynomial = polynomial
        self.table = self.generate_crc8_table()

    def generate_crc8_table(self):
        """Generate CRC8 table."""
        table = []
        for i in range(256):  # 8-bit values
            crc = i
            for j in range(8):  # 8-bit polynomial
                if crc & 0x80:  # if the MSB is set
                    crc = (crc << 1) ^ self.polynomial
                else:
                    crc <<= 1
                crc &= 0xFF
            table.append(crc)
        return table

    def compute(self, data):
        """Compute the CRC value of the input data."""
        crc = 0
        for byte in data:
            crc = self.table[crc ^ byte]
        return crc


class PacketGenerator:

    SFD_CODE = [
        0, 1, 0, 1, 0, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 1,
        1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 0, 0, 1, 0, 1, 0,
        0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0,
        0, 1, 1, 1, 1, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1
    ]

    PLCP_CODE = [
        1, 1, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 0,
        0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0,
        1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 0, 0, 1, 0, 0,
        1, 0, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0,
    ]

    S2P_MAPPING = {
        '0000': [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
        '0001': [1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0],
        '0010': [1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0],
        '0011': [1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1],
        '0100': [1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0],
        '0101': [1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1],
        '0110': [1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1],
        '0111': [1, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0],
        '1000': [1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0],
        '1001': [1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1, 0, 1],
        '1010': [1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1],
        '1011': [1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0],
        '1100': [1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1],
        '1101': [1, 0, 1, 0, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0],
        '1110': [1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0],
        '1111': [1, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1],

    }

    P2S_MAPPING = {tuple(v): k for k, v in S2P_MAPPING.items()}

    FSC_MAPPING = {
        0: [1, 0, 1, 0, 1, 0, 1, 0],
        1: [0, 1, 0, 1, 0, 1, 0, 1]
    }

    PILOT_INFO = {
        '000': -1,
        '001': -1,
        '010': 64,
        '011': 128,
        '100': -1,
        '101': -1,
        '110': np.inf,
        '111': -1
    }

    DATA_RATES = {
        '164kbps': '000',
        '328kbps': '001',
        '656kbps': '010',
        '13125kbps': '011',

    }
    DATA_RATE_MULTIPLIER = {
        '164kbps': 8,
        '328kbps': 4,
        '656kbps': 2,
        '13125kbps': 1,
    }

    SCRAMBLER_SEEDS = {
        0: "01101001010101000000000101010010",
        1: "10001010010111110110001000011111"
    }

    def __init__(self, spreading_enabled=False):
        self._spreading_enabled = spreading_enabled
        self.plcp_info = [0] * 32
        self.shift_register = [0] * 32
        self.scrambler_seed = 0
        self.set_scrambler_seed(self.scrambler_seed)
        self.PILOT_SEQUENCE = self.spread_with_fsc(self.SFD_CODE, 1)

    def set_scrambler_seed(self, seed):
        if seed in self.SCRAMBLER_SEEDS:
            self.shift_register = [int(i) for i in self.SCRAMBLER_SEEDS[seed]]
        else:
            raise ValueError(f"Invalid seed value: {seed}")

    def increment_scrambler_seed(self):
        self.scrambler_seed = 1 if self.scrambler_seed == 0 else 0
        self.set_scrambler_seed(self.scrambler_seed)

    def scrambler(self, input_data):
        scrambled_data = []
        for bit in input_data:
            # z[n] = z[n-11] mod2+ z[n-31] mod2+ z[n-32]
            feedback = (
                self.shift_register[21] ^ self.shift_register[1] ^ self.shift_register[0]) & 1
            output = feedback ^ int(bit)
            scrambled_data.append(output)
            self.shift_register.pop(0)
            self.shift_register.append(feedback)
        self.increment_scrambler_seed()

        return scrambled_data

    def symbol_to_chip(self, symbol):
        return self.S2P_MAPPING[symbol]

    def generate_PLCP_preamble(self):
        return self.spread_with_fsc(self.PLCP_CODE*4, 1)

    def generate_PLCP_header(self, data_rate, pilot_info, burst_mode, psdu_length):
        dt_rt = self.DATA_RATES.get(data_rate, '111')
        plcp_header = dt_rt + pilot_info + '00' + burst_mode + '00' + \
            str(self.scrambler_seed) + '0000' + format(psdu_length, '08b')
        plcp_header += header_crc_bits(plcp_header)
        return plcp_header  # Modify as needed

    def spread_with_fsc(self, orthogonal_code, data_rate):
        if not self._spreading_enabled:
            return orthogonal_code
        if data_rate == 0:
            return orthogonal_code
        spread_sequence = []
        for bit in orthogonal_code:
            spread_sequence.extend(
                SF_ZERO*data_rate if bit == 0 else SF_ONE*data_rate)
        return spread_sequence

    def data_to_chip(self, data_bits, data_rate):
        chip_data = []
        for i in range(0, len(data_bits), 4):
            symbol = ''.join(map(str, data_bits[i:i+4]))
            orthogonal_code = self.symbol_to_chip(symbol)
            chip_data.extend(self.spread_with_fsc(orthogonal_code, data_rate))

        return chip_data

    def generate_PSDU(self, data_bits, data_rate):
        return self.data_to_chip(data_bits, self.DATA_RATE_MULTIPLIER.get(data_rate, 0))

    def insert_pilot_sequence(self, data, pilot_info, data_rate):
        insertion_period = self.PILOT_INFO[pilot_info] * (self.DATA_RATE_MULTIPLIER.get(data_rate, 1)*SF_BIT_LEN - (
            self.DATA_RATE_MULTIPLIER.get(data_rate, 1)*SF_BIT_LEN-1) * (1-self._spreading_enabled)) * 32
        if insertion_period > len(data):
            return data
        if not insertion_period:
            return data

        # print(f'pilots in psdu{len(data)/(insertion_period )}')
        result = []
        for i in range(0, len(data), insertion_period):  # Convert octets to bits
            chunk = data[i:i + insertion_period]

            result.extend(chunk)
            # print(len(result))
            if(len(chunk)<insertion_period):
                break
            result.extend(self.PILOT_SEQUENCE)

        return result

    def scramble_PSDU(self, psdu):
        return self.scrambler(psdu)

    def generate_packet(self, data, data_rate, pilot_info, burst_mode, psdu_length):

        preamble = self.generate_PLCP_preamble()
        plcp_header = self.generate_PLCP_header(
            data_rate, pilot_info, burst_mode, psdu_length) #
        plcp_header = self.data_to_chip(
            plcp_header, self.DATA_RATE_MULTIPLIER.get('164kbps', 0))
        scrambled_psdu = self.scramble_PSDU(data)
        psdu_data = self.generate_PSDU(scrambled_psdu, data_rate)
        psdu_data_with_pilot = self.insert_pilot_sequence(
            psdu_data, pilot_info, data_rate)

        # Combine all to get the final packet
        return preamble + self.PILOT_SEQUENCE + plcp_header + psdu_data_with_pilot
    
    def generate_packet_no_preamble(self, data, data_rate, pilot_info, burst_mode, psdu_length):

        plcp_header = self.generate_PLCP_header(
            data_rate, pilot_info, burst_mode, psdu_length)
        plcp_header = self.data_to_chip(
            plcp_header, self.DATA_RATE_MULTIPLIER.get('164kbps', 0))
        scrambled_psdu = self.scramble_PSDU(data)
        psdu_data = self.generate_PSDU(scrambled_psdu, data_rate)
        psdu_data_with_pilot = self.insert_pilot_sequence(
            psdu_data, pilot_info, data_rate)

        # Combine all to get the final packet
        return plcp_header + psdu_data_with_pilot


class PacketParser:
    def __init__(self, packet_generator):
        self.packet_generator = packet_generator
    
    def detect_sequence(self, signal, code, threshold, DRAW=False):
            """Detects the code in the signal using cross-correlation and a threshold."""
            if len(signal) < len(code):
                return [0]
            if max(signal) >1 and min(signal) < 0:
                signal = (signal - min(signal))/(max(signal)-min(signal))
            if max(signal) <= 1 and min(signal) >= 0:
                signal = signal*2-1
            if max(code) <= 1 and min(code) >= 0:
                code = code*2-1
            #signal[signal == 0] = 1
            correlation_output = np.correlate(signal, code, mode='same')
            #print(signal[1000:1064], code)
            detections = np.where(correlation_output >= threshold)[0]
            if len(detections) == 0:
                detections = np.where(-correlation_output >= threshold)[0]
                if len(detections) == 0:
                    #DRAW = True
                    return []
            if DRAW:
                plt.figure(1)
                plt.plot(signal)
                plt.title('signal')
                plt.figure(2)
                plt.plot(correlation_output)
                plt.title('correlation_output')
                plt.show()
            return detections - len(code)//2  # Adjust indices to represent start of code
    
    def split_channel(self, signal, plcp_code, sdf_code, plcp_threshold, sdf_threshold):
            """Extract whole packets from a continuous data stream."""
            # Detect the preamble and SDF
            if len(signal) < 1000:
                return []
            preamble_starts = self.detect_sequence(signal, plcp_code, plcp_threshold, DRAW=False)
            sdf_starts = self.detect_sequence(signal, sdf_code, sdf_threshold)
            if len(preamble_starts) == 0 or len(sdf_starts) == 0:
                return [signal]
            delta = 64  
            groups = []
            current_group = [preamble_starts[0]]
            for i in range(1, len(preamble_starts)):
                # Check the delta
                if preamble_starts[i] - preamble_starts[i-1] > delta:
                    # Check if there is an element in array2 that is delta away from the last element of current_group
                    if current_group[-1] + delta in sdf_starts:
                        groups.append((current_group[-1] + delta*2,preamble_starts[i]))
                    current_group = []
                current_group.append(preamble_starts[i])


            if current_group[-1] + delta in sdf_starts:
                groups.append((current_group[0],len(signal)))

            return [np.array(signal[l:r]) for l,r in groups]
    
    def extract_packets(self, signal, plcp_code, sdf_code, plcp_threshold, sdf_threshold):
            """Extract whole packets from a continuous data stream."""
            # Detect the preamble and SDF
            if len(signal) < 1000:
                return []
            preamble_starts = self.detect_sequence(signal, plcp_code, plcp_threshold, DRAW=False)
            sdf_starts = self.detect_sequence(signal, sdf_code, sdf_threshold)
            if len(preamble_starts) == 0 or len(sdf_starts) == 0:
                return [signal]
            delta = 64  
            groups = []
            current_group = [preamble_starts[0]]
            for i in range(1, len(preamble_starts)):
                # Check the delta
                if preamble_starts[i] - preamble_starts[i-1] > delta:
                    # Check if there is an element in array2 that is delta away from the last element of current_group
                    if current_group[-1] + delta in sdf_starts:
                        groups.append((current_group[-1] + delta*2,preamble_starts[i]))
                    current_group = []
                current_group.append(preamble_starts[i])


            if current_group[-1] + delta in sdf_starts:
                groups.append((current_group[-1] + delta*2,len(signal)))

            return [np.array(signal[l:r]) for l,r in groups]


    def delete_preamble(self, channel):
        #preamble_start = np.argmax(np.correlate(channel, self.packet_generator.generate_PLCP_preamble()))
        preamble_start = 0
        return channel[preamble_start+len(self.packet_generator.generate_PLCP_preamble()):]

    def delete_sfd(self, channel):
        sfd_start = np.argmax(np.correlate(channel[:len(
            self.packet_generator.PILOT_SEQUENCE)*2], self.packet_generator.PILOT_SEQUENCE))

        return channel[sfd_start+len(self.packet_generator.PILOT_SEQUENCE):]

    def decode_PLCP_header(self, packet):
        # example plcp header - '000010000000000011111111'
        # example plcp header chips - '11111111111111111111111100000000111111111111111111111111111111111001011001101001100101100110100111111111000000001100001111000011'
        flipped = False
        
        if max(packet) == 1 and min(packet) == -1:
            packet = (packet+1)//2
        
        header_length = int(32/4*16) * (8*SF_BIT_LEN - ((8*SF_BIT_LEN-1)*(
            1-self.packet_generator._spreading_enabled)))  # 32 for header and 8 for CRC 8192
        plcp_header = packet[:header_length]
        
        if max(plcp_header) == 1 and min(plcp_header) == -1:
            plcp_header = (plcp_header+1)//2
        
        depspreaded_header = self.despread_data(
            plcp_header, self.packet_generator.DATA_RATE_MULTIPLIER.get('164kbps', 0)*SF_BIT_LEN)
        serialized_header = self.serialize_data(depspreaded_header)
        header = ''.join(serialized_header)
        
        crc_ok = header[-8:] == header_crc_bits(header[:-8])
        if not crc_ok:
            plcp_header = (1-plcp_header)
            depspreaded_header = self.despread_data(
                plcp_header, self.packet_generator.DATA_RATE_MULTIPLIER.get('164kbps', 0)*SF_BIT_LEN)
            serialized_header = self.serialize_data(depspreaded_header)
            header = ''.join(serialized_header)
            crc_ok = header[-8:] == header_crc_bits(header[:-8])
            if crc_ok:
                packet = (1-packet)
                
        # You can now decode each field of the header based on the structure you've defined
        # and return them.
        if len(header) < 4:
            return '164kbps', '110', '0', '0', 255, packet[header_length:], False
        data_rate = {code: rate for rate,
                     code in self.packet_generator.DATA_RATES.items()}.get(header[:3])
        if len(header) < 7:
            return data_rate, '110', '0', '0', 255, packet[header_length:], False
        pilot_info = header[3:6]
        if len(header) < 9:
            return data_rate, pilot_info, '0', '0', 255, packet[header_length:], False
        burst_mode = header[8]
        if len(header) < 12:
            return data_rate, pilot_info, burst_mode, '0', 255, packet[header_length:], False
        seed = header[11]
        if len(header) < 24:
            return data_rate, pilot_info, burst_mode, seed, 255, packet[header_length:], False
        psdu_length = int(header[16:-8], 2)
        if len(header) < 32:
            return data_rate, pilot_info, burst_mode, seed, psdu_length, packet[header_length:], False
        crc_ok = header[-8:] == header_crc_bits(header[:-8])

        return data_rate, pilot_info, burst_mode, seed, psdu_length, packet[header_length:], crc_ok

    def remove_pilot_sequences(self, packet, pilot_info, pilot_count):

        def remove_pilots(payload, sdf_code, sdf_threshold):
            """Remove pilot sequences from the payload."""
            # Detect the pilot sequences
            if len(payload) < len(sdf_code):
                return payload
            pilot_starts = self.detect_sequence(payload, sdf_code, sdf_threshold)

            cleaned_payload = np.copy(payload)

            # Remove detected pilots
            indices_to_delete = [np.arange(pilot_start, pilot_start + len(sdf_code)) for pilot_start in pilot_starts]
            cleaned_payload = np.delete(cleaned_payload, indices_to_delete)

            return cleaned_payload
        return remove_pilots(np.array(packet), np.array(self.packet_generator.PILOT_SEQUENCE), 50)


    # [0,1,] => 1
    # [0,1,0,1,] => 1
    # [0,1,0,1,0,1,0,1,] => 1

    def despred(self, byte_seq):
        return np.sign(np.dot(byte_seq*2-1, np.array([-1, 1]*int(len(byte_seq)/2))))
        # return all(np.array([
        #        np.rint(statistics.median(byte_seq[::2])+np.finfo(np.float64).resolution),
        #        np.rint(statistics.median(byte_seq[1::2])+np.finfo(np.float64).resolution)
        #        ]) == np.array([0,1]))

    def despread_data(self, packet, data_rate):
        # The reverse of spread_with_fsc and symbol_to_chip
        if not self.packet_generator._spreading_enabled:
            return packet
        reshaped_packet = packet.reshape(
            (int(len(packet)/(data_rate)), (data_rate)))

        return np.fromiter(map(self.despred, reshaped_packet), dtype=int)
        # despread_sequence = np.empty(len(reshaped_packet),dtype=int)
        # for indx,symbol in enumerate(reshaped_packet):
        #    #symbol = np.array([int(bit) for bit in packet[i:i+2*data_rate]])
        #    if len(symbol) < 2*data_rate:
        #        break
        #    sym_temp_val = np.array([
        #        int(np.rint(statistics.median(symbol[::2])+np.finfo(np.float64).resolution)),
        #        int(np.rint(statistics.median(symbol[1::2])+np.finfo(np.float64).resolution))
        #        ])
        #    sym_sum =  sum(sym_temp_val)
        #
        #    if sym_sum!= 1:
        #        symbol = int(sum(np.abs(symbol - np.array([0,1]*data_rate))) < sum(np.abs(symbol - np.array([1,0]*data_rate))))
        #    else:
        #        symbol = int(all(sym_temp_val == np.array([0,1])))
        #    despread_sequence[indx] = symbol
        # return despread_sequence

    def serialize_data(self, packet):
       # Define the size of the Walsh matrix
        data_bits = []
        for i in range(0, len(packet), 16):
            # math black magic with walsh matrix
            try:
                symbol = '{0:04b}'.format(
                    np.argmax(np.dot((np.array(packet[i:i+16])-0.5)*2, np.array(hadamard_matrix))))
            except:
                symbol = '0000'

            # symbol = self.packet_generator.P2S_MAPPING[tuple(packet[i:i+16])]
            data_bits.extend(symbol)
        return data_bits
        # return None

    def descramble(self, scrambled_data, seed):
        # Setting seed for descrambling
        self.packet_generator.set_scrambler_seed(int(seed))
        return self.packet_generator.scrambler(scrambled_data)

    def unpack_packet(self, channel):
        # Unpacking based on the mentioned structure

        #packet = self.delete_preamble(channel)
        #packet = self.delete_sfd(packet)
        data_rate, pilot_info, burst_mode, seed, psdu_length, remaining_packet, crc_ok = self.decode_PLCP_header(
            channel)
        # print(f'data rate = {data_rate}, pilot info = {pilot_info}, burst mode = {burst_mode}, seed = {seed}, psdu length = {psdu_length}, crc ok = {crc_ok}')
        if not crc_ok:
            print('CRC failed')
        packet_length_without_pilot = psdu_length * 8 / 4 * 16 * (self.packet_generator.DATA_RATE_MULTIPLIER.get(data_rate, 1)*SF_BIT_LEN - (
            self.packet_generator.DATA_RATE_MULTIPLIER.get(data_rate, 1)*SF_BIT_LEN-1) * (1-self.packet_generator._spreading_enabled))
        insertion_period = self.packet_generator.PILOT_INFO[pilot_info] * (self.packet_generator.DATA_RATE_MULTIPLIER.get(data_rate, 1)*SF_BIT_LEN - (
            self.packet_generator.DATA_RATE_MULTIPLIER.get(data_rate, 1)*SF_BIT_LEN-1) * (1-self.packet_generator._spreading_enabled)) * 32
        # pilot_count = int(packet_length_without_pilot/(self.packet_generator.PILOT_INFO[pilot_info]*8*8))
        pilot_count = int(packet_length_without_pilot//insertion_period)
        packet_length_with_pilot = int(
            pilot_count*len(self.packet_generator.PILOT_SEQUENCE) + packet_length_without_pilot)
        data_with_pilots = remaining_packet[:packet_length_with_pilot]
        data_without_pilots = self.remove_pilot_sequences(
            data_with_pilots, pilot_info, pilot_count)
        despreaded_data = self.despread_data(
            data_without_pilots, self.packet_generator.DATA_RATE_MULTIPLIER.get(data_rate, 1)*SF_BIT_LEN)
        serialized_data = self.serialize_data(despreaded_data)
        descrambled_data = self.descramble(serialized_data, seed)

        return ''.join(map(str, descrambled_data)),data_rate, pilot_info, burst_mode, seed, psdu_length,crc_ok,data_without_pilots

    def bin_to_str(self, x):
        ''' Converts a Binary String to an (ASCII) string'''
        my_str = [chr(int(x[i:i+8], 2)) for i in range(0, len(x), 8)]
        my_str = ''.join(my_str)
        return my_str


# # # Example usage

def find_best_match(received_text, original_parts):
    best_ratio = 0
    best_match = r'\0'
    for part in original_parts:
        # Use difflib to find the closest match to the received_text
        s = difflib.SequenceMatcher(None, received_text, part)
        if s.ratio() > best_ratio:
            best_ratio = s.ratio()
            best_match = part
    return best_match

if __name__ == "__main__":

    import numpy as np

# Example 64-bit normalized Gold code
    


    generator = PacketGenerator(spreading_enabled=False)
    raw_data_to_send = 'We are going to die, and that makes us the lucky ones. Most people are never going to die because they are never going to be born. The potential people who could have been here in my place but who will in fact never see the light of day outnumber the sand grains of Arabia. Certainly those unborn ghosts include greater poets than Keats, scientists greater than Newton. We know this because the set of possible people allowed by our DNA so massively exceeds the set of actual people. In the teeth of these stupefying odds it is you and I, in our ordinariness, that are here. We privileged few, who won the lottery of birth against all odds, how dare we whine at our inevitable return to that prior state from which the vast majority have never stirred?'
    multiple_data_to_send = [raw_data_to_send[i:i+255]
                             for i in range(0, len(raw_data_to_send), 255)]
    channel = []
    for indx, data_to_send in enumerate(multiple_data_to_send):
        my_data = [format(ord(i), '08b') for i in data_to_send]
        my_data = ''.join(my_data)
        packet = np.array(generator.generate_packet(data=my_data, data_rate='13125kbps',
                          pilot_info='010', burst_mode='0', psdu_length=len(my_data)//8))
        empty_channel = np.array([0] * 1000)
        # print(f'data = {my_data}')
        # print(f'packet {indx} length = {len(packet)}')
        # print(''.join([str(i) for i in packet]))
        channel = np.concatenate(
            [channel, empty_channel, packet, empty_channel])

    print('Sent data: ', raw_data_to_send)

    with open('data_to_send', 'wb') as file:
        np.save(file, channel.astype(np.uint8))

    parser = PacketParser(generator)
    channel = channel.astype(int)
    pckts = parser.count_packets(channel)
    
    data = []

    #for pos_l, pos_r in zip(positions[:-1], positions[1:]):
    for pckt in pckts:
        unpacked_data = parser.unpack_packet(pckt)

        print('Received data: ', parser.bin_to_str(unpacked_data))
